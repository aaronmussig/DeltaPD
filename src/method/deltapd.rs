use crate::method::common_taxa::{find_common_ids_in_pdms, CommonTaxa};
use crate::model::error::{DeltaPDError, DeltaPDResult};
use crate::model::job::{Job, JobData};
use crate::model::linalg::{LinearModel, LinearModelNew};
use crate::model::metadata::MetadataFile;
use crate::model::output::DeltaPdOutput;
use crate::model::params::Params;
use crate::model::pdm::{QryDistMatrix, RefDistMatrix};
use crate::model::result::OutputResultSmall;
use crate::stats::jackknife::{calc_prop_std_error_of_jackknife, calc_relative_jackknife_influence};
use indicatif::ParallelProgressIterator;
use phylodm::tree::Taxon;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::prelude::*;
use std::time::Instant;

use crate::stats::vec::{calc_mean, calc_median, calc_stddev};

/// The main method that runs DeltaPD
pub fn run_deltapd(qry_mat: &QryDistMatrix, ref_mat: &RefDistMatrix, metadata_file: &MetadataFile, params: &Params) -> DeltaPDResult<Vec<OutputResultSmall>> {

    // 1. Identify the common taxa in both trees
    let common_ids = find_common_ids_in_pdms(&ref_mat.dm.taxa, &qry_mat.dm.taxa, &metadata_file)?;

    // 2. Generate the query taxa that will be processed
    let qry_taxa_to_process = generate_taxa_to_process(params, &common_ids)?;

    // 3. Generate a queue of jobs to be processed
    let jobs = create_jobs(&qry_taxa_to_process, params)?;
    println!("Created {} jobs.", jobs.len());

    // 4. Process each job
    let previous_time = Instant::now();
    let dpd_output: Vec<DeltaPdOutput> = if params.cpus > 1 {
        println!("Running with {} cpus", params.cpus);
        jobs.par_iter().progress_count(jobs.len() as u64).map(|cur_job| {
            run_deltapd_on_taxon_test(cur_job, qry_mat, ref_mat, metadata_file, params).unwrap()
        }).flatten().collect()
    } else {
        println!("Running with 1 CPU.");
        jobs.iter().map(|cur_job| {
            run_deltapd_on_taxon_test(cur_job, qry_mat, ref_mat, metadata_file, params).unwrap()
        }).flatten().collect()
    };
    println!("Ran main method in {:?}", previous_time.elapsed());

    // 5. Transform the results into the expected output format
    let results = transform_results_to_output(&dpd_output, qry_taxa_to_process.len())?;

    // 7. Return the results
    Ok(results)
}


/// Given the parameters and the common taxa, generate the taxa that will be processed
/// If no taxa are provided in the params, then all common taxa will be processed.
pub fn generate_taxa_to_process<'a>(params: &Params, common_ids: &'a CommonTaxa) -> DeltaPDResult<HashSet<&'a Taxon>> {
    if params.taxa.is_empty() {
        if common_ids.qry_taxa.len() == 0 {
            Err(DeltaPDError::Error("No query taxa mapped to their reference taxon.".to_string()))
        } else {
            Ok(common_ids.qry_taxa.clone())
        }
    } else {
        let mut out = HashSet::with_capacity(params.taxa.len());
        for taxon in &params.taxa {
            let taxon_obj = Taxon(taxon.to_string());
            let taxon_ref = common_ids.qry_taxa.get(&taxon_obj).unwrap();
            out.insert(*taxon_ref);
        }
        if out.is_empty() {
            Err(DeltaPDError::Error("No query taxa mapped to their reference taxon.".to_string()))
        } else {
            Ok(out)
        }
    }
}


/// Create a list of taxa to be processed and their replicates.
pub fn create_jobs<'a>(qry_taxa: &'a HashSet<&Taxon>, params: &Params) -> DeltaPDResult<Vec<Job<'a>>> {
    let n_jobs = params.replicates * qry_taxa.len();
    let mut out = Vec::with_capacity(n_jobs);
    for &taxon in qry_taxa {
        for replicate_id in 0..params.replicates {
            out.push(
                Job::new(
                    taxon,
                    replicate_id,
                )
            )
        }
    }
    if out.is_empty() {
        Err(DeltaPDError::Error("No jobs created.".to_string()))
    } else {
        Ok(out)
    }
}

/// The main method that processes a single job.
pub fn run_deltapd_on_taxon_test(cur_job: &Job, qry_mat: &QryDistMatrix, ref_mat: &RefDistMatrix, metadata_file: &MetadataFile, params: &Params) -> DeltaPDResult<Vec<DeltaPdOutput>> {

    // Extract the query taxa of interest
    let q_taxon = cur_job.taxon;

    // Get the X/Y data for the current job based on subsampling
    let job_data = JobData::new(cur_job, qry_mat, ref_mat, metadata_file, params);
    let ref_data = job_data.ref_data.as_slice();
    let qry_data = job_data.qry_data.as_slice();

    let linmodel = LinearModelNew::new(
        ref_data,
        qry_data,
    )?;

    // Calculate the base model
    // let base_model = linmodel.compute(None);
    // println!("{:?} Base model: {:?}", q_taxon, base_model);

    // Calculate it again using the main method
    // let true_model = calc_repeated_median(&ref_data, &qry_data);
    // // println!("{:?} True model gradient: {:?}, intercept: {:?}", q_taxon,
    //          true_model.gradient == base_model.gradient,
    //             true_model.intercept == base_model.intercept
    // );

    // Do the jackknifing and test that compute works
    let mut model_eval: Vec<(&Taxon, LinearModel)> = Vec::with_capacity(job_data.label_bitmask.len());
    for (cur_taxon, cur_taxon_bv) in job_data.label_bitmask {
        let cur_bitvec = Some(&cur_taxon_bv);
        let cur_model_params = linmodel.compute(cur_bitvec);
        let cur_model_eval = linmodel.fit(cur_bitvec, &cur_model_params, params.model_error, params.model_corr);
        let model_output = LinearModel {
            params: cur_model_params,
            eval: cur_model_eval,
        };
        model_eval.push((cur_taxon, model_output));
    }

    // Calculate the relative influence function
    let relative_influence = calc_relative_jackknife_influence(
        &model_eval
            .iter()
            .map(|x| x.1.eval.error)
            .collect::<Vec<f64>>(),
    );


    let std_error = calc_prop_std_error_of_jackknife(&relative_influence);

    // Collect the results
    let mut output: Vec<DeltaPdOutput> = Vec::with_capacity(model_eval.len());
    for (i, (cur_taxon, cur_model)) in model_eval.iter().enumerate() {
        let cur_rinf = relative_influence[i];
        let cur_err = std_error[i];
        let cur_model_grad = cur_model.params.gradient;
        let cur_model_int = cur_model.params.intercept;

        // Create the output object for those we are interested in
        if params.taxa.is_empty() {
            let cur_output = DeltaPdOutput::new(
                cur_taxon,
                cur_job.replicate,
                cur_rinf,
                cur_err,
                cur_model_grad,
                cur_model_int,
            );
            output.push(cur_output);
        } else {
            if params.taxa.contains(&cur_taxon.0) {
                let cur_output = DeltaPdOutput::new(
                    cur_taxon,
                    cur_job.replicate,
                    cur_rinf,
                    cur_err,
                    cur_model_grad,
                    cur_model_int,
                );
                output.push(cur_output);
            }
        }
    }


    // If we're in debug mode then output the X/Y coordinates from, this taxon to disk
    if params.debug {
        {
            let mut file_path = params.output_dir.clone();
            file_path.push(format!("xy_{}_{}.tsv", q_taxon.0, cur_job.replicate));

            let mut file = File::create(file_path).unwrap();
            file.write_all(b"query_taxon_i\tqry_taxon_j\tref_taxon_i\tqry_taxon_j\tquery_dist\tref_dist\n").unwrap();
            for i in 0..qry_data.len() {
                let line = format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\n",
                    job_data.qry_labels[i].0.0,
                    job_data.qry_labels[i].1.0,
                    job_data.ref_labels[i].0.0,
                    job_data.ref_labels[i].1.0,
                    qry_data[i],
                    ref_data[i]
                );
                file.write_all(line.as_bytes()).unwrap();
            }
        }
        {
            let mut file_path = params.output_dir.clone();
            file_path.push(format!("model_{}_{}.tsv", q_taxon.0, cur_job.replicate));

            let mut file = File::create(file_path).unwrap();
            file.write_all(b"taxon\trelative_influence\tstd_error\tgradient\tintercept\terr\tcorr\n").unwrap();
            for (i, (cur_taxon, cur_model)) in model_eval.iter().enumerate() {
                let cur_rinf = relative_influence[i];
                let cur_err = std_error[i];
                let cur_model_grad = cur_model.params.gradient;
                let cur_model_int = cur_model.params.intercept;
                let cur_error = cur_model.eval.error;
                let cur_corr = cur_model.eval.corr;
                let line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\n", cur_taxon.0, cur_rinf, cur_err, cur_model_grad, cur_model_int, cur_error, cur_corr);
                file.write_all(line.as_bytes()).unwrap();
            }
        }
    }

    // Done, return the output
    Ok(output)
}


pub fn transform_results_to_output(dpd_output: &[DeltaPdOutput], n_query_taxa: usize) -> DeltaPDResult<Vec<OutputResultSmall>> {
    let mut out: Vec<OutputResultSmall> = Vec::with_capacity(n_query_taxa);

    // Create a hashmap to store the values of interest
    let mut taxon_to_data: HashMap<&Taxon, Vec<f64>> = HashMap::with_capacity(n_query_taxa);
    for cur_job in dpd_output {
        let cur_taxon = &cur_job.taxon;
        let cur_data = taxon_to_data.entry(cur_taxon).or_insert_with(|| Vec::with_capacity(cur_job.replicate + 1));
        cur_data.push(cur_job.prop_std_error);
    }

    // Transform the results into the expected output format
    for (cur_taxon, cur_data) in taxon_to_data {

        out.push(OutputResultSmall {
            taxon: cur_taxon.clone(),
            error_mean: calc_mean(&cur_data),
            error_median: calc_median(&cur_data),
            error_std: calc_stddev(&cur_data),
        });
    }

    Ok(out)
}

