use crate::method::common_taxa::{find_common_ids_in_pdms, CommonTaxa};
use crate::method::knn::create_vecs_from_knn3;
use crate::model::error::{DeltaPDError, DeltaPDResult};
use crate::model::job::{Job, JobData};
use crate::model::linalg::{LinearModel, LinearModelNew};
use crate::model::metadata::MetadataFile;
use crate::model::params::Params;
use crate::model::pdm::{QryDistMatrix, RefDistMatrix};
use crate::model::result::OutputResultSmall;
use crate::stats::jackknife::{calc_prop_std_error_of_jackknife, calc_relative_jackknife_influence, jackknife_fn};
use crate::stats::linalg::calc_repeated_median;
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use phylodm::tree::Taxon;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use rayon::prelude::*;
use std::collections::HashSet;
use std::fs::File;
use std::io::prelude::*;
use std::time::Instant;

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
            let taxon_obj_2 = common_ids.qry_taxa.get(&taxon_obj).unwrap();
            out.insert(*taxon_obj_2);
        }
        if out.is_empty() {
            Err(DeltaPDError::Error("No query taxa mapped to their reference taxon.".to_string()))
        } else {
            Ok(out)
        }
    }
}

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

pub fn run_deltapd(qry_mat: &QryDistMatrix, ref_mat: &RefDistMatrix, metadata_file: &MetadataFile, params: &Params) -> DeltaPDResult<Vec<OutputResultSmall>> {

    // Identify the common taxa between both trees
    let previous_time = Instant::now();
    let common_ids = find_common_ids_in_pdms(&ref_mat.dm.taxa, &qry_mat.dm.taxa, &metadata_file)?;
    println!("Found common IDs in {:?}", previous_time.elapsed());

    // Generate the list of query taxa that should be processed, all if no filtering is required
    let qry_taxa_to_process = generate_taxa_to_process(params, &common_ids)?;

    // Create a queue of jobs for processing
    let previous_time = Instant::now();
    let jobs = create_jobs(&qry_taxa_to_process, params)?;
    println!("Created {:?} jobs in {:?}", jobs.len(), previous_time.elapsed());

    println!("Iterating over each taxon (par)");
    let previous_time = Instant::now();

    if params.cpus > 1 {
        println!("Running with {} cpus", params.cpus);
        let x: Vec<_> = jobs.par_iter().progress_count(jobs.len() as u64).map(|(cur_job)| {
            run_deltapd_on_taxon_test(cur_job, qry_mat, ref_mat, metadata_file, params);
        }).collect();
        // let x: Vec<_> = create_pool(params.cpus)?.install(|| {
        //     jobs.iter().map(|(cur_job)| {
        //         run_deltapd_on_taxon_test(cur_job, qry_mat, ref_mat, metadata_file, params);
        //     })
        // }).collect();
    } else {
        let x: Vec<_> = jobs.iter().map(|(cur_job)| {
            run_deltapd_on_taxon_test(cur_job, qry_mat, ref_mat, metadata_file, params);
        }).collect();
    }
    println!("Ran main method (par) in {:?}", previous_time.elapsed());


    panic!("Stop!");


    println!("Iterating over each taxon (par)");
    let previous_time = Instant::now();
    let output: Vec<OutputResultSmall> = jobs.par_iter().map(|(cur_job)| {
        run_deltapd_on_taxon(cur_job, qry_mat, ref_mat, metadata_file, params).unwrap()
    }).collect();
    println!("Ran main method (par) in {:?}", previous_time.elapsed());

    Ok(output)
}

pub fn run_deltapd_on_taxon_test(cur_job: &Job, qry_mat: &QryDistMatrix, ref_mat: &RefDistMatrix, metadata_file: &MetadataFile, params: &Params) {
    // Extract the query taxa of interest
    let q_taxon = cur_job.taxon;

    // Get the X/Y data for the current job based on subsampling
    let job_data = JobData::new(cur_job, qry_mat, ref_mat, metadata_file, params.sample_size);
    let ref_data = job_data.ref_data.as_slice();
    let qry_data = job_data.qry_data.as_slice();

    let linmodel = LinearModelNew::new(
        ref_data,
        qry_data,
    );

    // Calculate the base model
    let base_model = linmodel.compute(None);
    // println!("{:?} Base model: {:?}", q_taxon, base_model);

    // Calculate it again using the main method
    let true_model = calc_repeated_median(&ref_data, &qry_data);
    // // println!("{:?} True model gradient: {:?}, intercept: {:?}", q_taxon,
    //          true_model.gradient == base_model.gradient,
    //             true_model.intercept == base_model.intercept
    // );

    // Do the jackknifing and test that compute works
    let mut model_eval: Vec<(&Taxon, LinearModel)> = Vec::with_capacity(job_data.label_bitmask.len());
    for (cur_taxon, cur_taxon_bv) in job_data.label_bitmask {
        let cur_bitvec = Some(&cur_taxon_bv);
        let cur_model_params = linmodel.compute(cur_bitvec);
        let cur_model_eval = linmodel.fit(cur_bitvec, cur_model_params.clone(), params.model_error, params.model_corr);
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
    let std_error_total = std_error.iter().sum::<f64>();

    // As I am lazy write the data to disk
    let mut file = File::create(format!("/tmp/dpd/deltapd_xy_{}_{}.tsv", q_taxon.0, cur_job.replicate)).unwrap();
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

    // Write out the relative influence values
    let mut file = File::create(format!("/tmp/dpd/deltapd_relative_influence_{}_{}.tsv", q_taxon.0, cur_job.replicate)).unwrap();
    file.write_all(b"taxon\trelative_influence\tstd_error\tgradient\tintercept\n").unwrap();
    for (i, (cur_taxon, cur_model)) in model_eval.iter().enumerate() {
        let cur_rinf = relative_influence[i];
        let cur_err = std_error[i];
        let cur_model_grad = cur_model.params.gradient;
        let cur_model_int = cur_model.params.intercept;
        let line = format!("{}\t{}\t{}\t{}\t{}\n", cur_taxon.0, cur_rinf, cur_err, cur_model_grad, cur_model_int);
        file.write_all(line.as_bytes()).unwrap();
    }
}


pub fn run_deltapd_on_taxon(cur_job: &Job, qry_mat: &QryDistMatrix, ref_mat: &RefDistMatrix, metadata_file: &MetadataFile, params: &Params) -> DeltaPDResult<OutputResultSmall> {
    let q_taxon = cur_job.taxon;

    // Create the XY vectors containing the data points from the KNN set
    let knn_as_linalg = create_vecs_from_knn3(&qry_mat, &ref_mat, q_taxon, metadata_file, params)?;

    // // Fit a model without any jackknifing
    // let linear_model_base = LinearModel::fit(
    //     params.model,
    //     params.model_error,
    //     params.model_corr,
    //     &knn_as_linalg.ref_data,
    //     &knn_as_linalg.qry_data,
    // );

    // Perform jacknifing on the data
    let jackknife_models = jackknife_fn(
        &knn_as_linalg.ref_data,
        &knn_as_linalg.qry_data,
        params.model,
        params.model_error,
        params.model_corr,
    )?;

    // Calculate the relative influence function
    let relative_influence = calc_relative_jackknife_influence(
        &jackknife_models
            .iter()
            .map(|x| x.eval.error)
            .collect::<Vec<f64>>(),
    );

    let std_error = calc_prop_std_error_of_jackknife(&relative_influence);

    Ok(OutputResultSmall {
        taxa: knn_as_linalg.qry_labels,
        std_error,
    })
    // Ok(OutputResult {
    //     linear_model_base,
    //     linear_model_jk: jackknife_models,
    //     relative_influence,
    //     std_error,
    //     query_taxon: q_taxon.clone(),
    //     knn_as_linalg,
    // })

}