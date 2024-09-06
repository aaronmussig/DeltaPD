use crate::method::common_taxa::{find_common_ids_in_pdms, CommonTaxa};
use crate::method::knn::{create_vecs_from_knn3, sample_with_replacement, sample_with_replacement_include};
use crate::model::error::DeltaPDResult;
use crate::model::job::Job;
use crate::model::metadata::MetadataFile;
use crate::model::params::Params;
use crate::model::pdm::{QryDistMatrix, RefDistMatrix};
use crate::model::result::OutputResultSmall;
use crate::stats::jackknife::{calc_prop_std_error_of_jackknife, calc_relative_jackknife_influence, jackknife_fn, jackknife_fn_masked};
use phylodm::tree::Taxon;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::time::Instant;
use std::fs::File;
use std::io::prelude::*;
use bitvec::bitvec;
use bitvec::vec::BitVec;
use crate::model::linalg::{LinearModel, LinearModelNew};
use indicatif::ParallelProgressIterator;
use rayon::iter::{ParallelIterator, IntoParallelRefIterator};
use crate::stats::linalg::{calc_repeated_median, RepeatedMedian};

pub fn generate_taxa_to_process<'a>(params: &Params, common_ids: &'a CommonTaxa) -> HashSet<&'a Taxon> {
    if params.taxa.is_empty() {
        common_ids.qry_taxa.clone()
    } else {
        let mut out = HashSet::with_capacity(params.taxa.len());
        for taxon in &params.taxa {
            let taxon_obj = Taxon(taxon.to_string());
            let taxon_obj_2 = common_ids.qry_taxa.get(&taxon_obj).unwrap();
            out.insert(*taxon_obj_2);
        }
        out
    }
}

pub fn create_jobs<'a>(qry_taxa: &'a HashSet<&Taxon>, params: &Params) -> Vec<Job<'a>> {
    let n_jobs = params.replicates * qry_taxa.len();
    let mut out = Vec::with_capacity(n_jobs);
    for taxon in qry_taxa {
        for replicate_id in 0..params.replicates {
            out.push(
                Job::new(
                    *taxon,
                    replicate_id,
                )
            )
        }
    }
    out
}

pub fn run_deltapd(qry_mat: &QryDistMatrix, ref_mat: &RefDistMatrix, metadata_file: &MetadataFile, params: &Params) -> DeltaPDResult<Vec<OutputResultSmall>> {

    // Identify the common taxa between both trees
    let previous_time = Instant::now();
    let common_ids = find_common_ids_in_pdms(&ref_mat.dm.taxa, &qry_mat.dm.taxa, &metadata_file)?;
    println!("Found common IDs in {:?}", previous_time.elapsed());

    // Generate the list of query taxa that should be processed, all if no filtering is required
    let qry_taxa_to_process = generate_taxa_to_process(params, &common_ids);

    // Create a queue of jobs for processing
    let previous_time = Instant::now();
    let jobs = create_jobs(&qry_taxa_to_process, params);
    println!("Created {:?} jobs in {:?}", jobs.len(), previous_time.elapsed());

    println!("Iterating over each taxon (par)");
    let previous_time = Instant::now();
    let x: Vec<_> = jobs.iter().map(|(cur_job)| {
        run_deltapd_on_taxon_test(cur_job, qry_mat, ref_mat, metadata_file, params);
    }).collect();
    // let x: Vec<_> = jobs.par_iter().progress_count(jobs.len() as u64).map(|(cur_job)| {
    //     run_deltapd_on_taxon_test(cur_job, qry_mat, ref_mat, metadata_file, params);
    // }).collect();
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
    /*
    The goal for this is to pick a specific taxon, then get the x nearest neighbours
    for the x nearest neighbours, get all the pairwise distances for these
    and use them to create the model

    from that we want to then calculate the relative influence of removing that specific taxon
     */

    let q_taxon = cur_job.taxon;

    // Sample taxa with replacement
    let qry_taxa_sampled = qry_mat.dm.sample_taxa_with_replacement_include(params.sample_size, q_taxon);
    let qry_taxa_sampled_unq = qry_taxa_sampled.iter().cloned().collect::<HashSet<&Taxon>>();

    let query_taxon_idx = *(qry_mat.dm.taxon_to_idx.get(q_taxon).unwrap());
    let query_ids_to_use = sample_with_replacement_include(qry_mat.dm.taxa.len(), params.sample_size, query_taxon_idx);
    let query_ids_to_use_unq = query_ids_to_use.iter().cloned().collect::<HashSet<usize>>();

    // Calculate the number of data points that will be included in this
    let n_data_points = query_ids_to_use.len() * (query_ids_to_use.len() - 1) / 2;

    // Collect the pairwise distances for those query taxa
    let mut qry_data: Vec<f64> = Vec::with_capacity(n_data_points);
    let mut ref_data: Vec<f64> = Vec::with_capacity(n_data_points);

    let mut qry_labels_i: Vec<String> = Vec::with_capacity(n_data_points);
    let mut qry_labels_j: Vec<String> = Vec::with_capacity(n_data_points);
    let mut ref_labels_i: Vec<String> = Vec::with_capacity(n_data_points);
    let mut ref_labels_j: Vec<String> = Vec::with_capacity(n_data_points);
    let mut mask: Vec<(usize, usize)> = Vec::with_capacity(n_data_points);
    let mut label_bitmask: HashMap<&Taxon, BitVec> = HashMap::with_capacity(query_ids_to_use_unq.len());

    // Initialise the bitmask for each taxon
    for qry_idx in query_ids_to_use_unq {
        let qry_taxon = &qry_mat.dm.taxa[qry_idx];
        label_bitmask.insert(qry_taxon, bitvec![0; n_data_points]);
    }

    // Iterate over each randomly sampled taxon index
    for (idx_i, &qry_idx_i) in query_ids_to_use.iter().enumerate() {

        // Collect the query and reference taxa (from)
        let qry_taxon_i = &qry_mat.dm.taxa[qry_idx_i];
        let ref_taxon_i_str = metadata_file.get_ref_taxon(&qry_taxon_i.0).unwrap();
        let ref_taxon_i = ref_mat.dm.get_taxon_from_string(ref_taxon_i_str);

        // Iterate over the lower triangle of the randomly sampled taxa indices
        for idx_j in 0..idx_i {

            // Collect the query and reference taxa (to)
            let qry_idx_j = query_ids_to_use[idx_j];
            let qry_taxon_j = &qry_mat.dm.taxa[qry_idx_j];
            let ref_taxon_j_str = metadata_file.get_ref_taxon(&qry_taxon_j.0).unwrap();
            let ref_taxon_j = ref_mat.dm.get_taxon_from_string(ref_taxon_j_str);

            // Get the distance between the two taxa
            let qry_dist = qry_mat.dm.distance(qry_taxon_i, qry_taxon_j);
            let ref_dist = ref_mat.dm.distance(ref_taxon_i, ref_taxon_j);

            // Store the values
            mask.push((qry_idx_i, qry_idx_j));
            label_bitmask.get_mut(qry_taxon_i).unwrap().set(idx_j, true);
            label_bitmask.get_mut(qry_taxon_j).unwrap().set(idx_i, true);

            qry_data.push(qry_dist);
            qry_labels_i.push(qry_taxon_i.0.clone());
            qry_labels_j.push(qry_taxon_j.0.clone());

            ref_data.push(ref_dist);
            ref_labels_i.push(ref_taxon_i.0.clone());
            ref_labels_j.push(ref_taxon_j.0.clone());
        }
    }

    let previous_time = Instant::now();
    let linmodel = LinearModelNew::new(
        &ref_data,
        &qry_data,
        params.model_error,
        params.model_corr
    );
    println!(">>>> Linear model new init {:?}", previous_time.elapsed());

    let previous_time = Instant::now();
    let cur_mask = label_bitmask.get(q_taxon).unwrap();
    let res1 = linmodel.compute(None);
    println!(">>>> Linear model new compute {:?}", previous_time.elapsed());

    let previous_time = Instant::now();
    let repeated_median = RepeatedMedian::new(&ref_data, &qry_data);
    let res2 = repeated_median.compute(&HashSet::new());
    println!(">>>> Repeated median class took {:?}", previous_time.elapsed());

    let previous_time = Instant::now();
    let res3 = repeated_median.compute(&HashSet::new());
    println!(">>>> Repeated median class (again) took {:?}", previous_time.elapsed());

    let previous_time = Instant::now();
    let res4 = calc_repeated_median(&ref_data, &qry_data);
    println!(">>>> Repeated median function took {:?}", previous_time.elapsed());

    let linear_model_base = LinearModel::fit(
        params.model,
        params.model_error,
        params.model_corr,
        &ref_data,
        &qry_data,
    );

    // Perform jacknifing on the data
    // When doing this do I want to remove ALL points that belong to a given taxon?
    let (jackknife_models, taxa_order) = jackknife_fn_masked(
        params.model,
        params.model_error,
        params.model_corr,
        &mask,
        &repeated_median
    ).unwrap();
    //
    // // As I am lazy write the data to disk
    // let mut file = File::create(format!("/tmp/dpd/deltapd_xy_{}_{}.tsv", q_taxon.0, cur_job.replicate)).unwrap();
    // file.write_all(b"query_taxon_i\tqry_taxon_j\tref_taxon_i\tqry_taxon_j\tquery_dist\tref_dist\n").unwrap();
    // for i in 0..qry_data.len() {
    //     let line = format!("{}\t{}\t{}\t{}\t{}\t{}\n", qry_labels_i[i], qry_labels_j[i], ref_labels_i[i], ref_labels_j[i], qry_data[i], ref_data[i]);
    //     file.write_all(line.as_bytes()).unwrap();
    // }
    //
    // // Calculate the relative influence function
    // let relative_influence = calc_relative_jackknife_influence(
    //     &jackknife_models
    //         .iter()
    //         .map(|x| x.eval.error)
    //         .collect::<Vec<f64>>(),
    // );
    //
    // let std_error = calc_prop_std_error_of_jackknife(&relative_influence);
    //
    // // Write out the relative influence values
    // let mut file = File::create(format!("/tmp/dpd/deltapd_relative_influence_{}_{}.tsv", q_taxon.0, cur_job.replicate)).unwrap();
    // file.write_all(b"taxon\trelative_influence\tstd_error\n").unwrap();
    // for (i, cur_taxon_idx) in taxa_order.iter().enumerate() {
    //     let cur_taxon = &qry_mat.dm.taxa[*cur_taxon_idx];
    //     let cur_rinf = relative_influence[i];
    //     let cur_err = std_error[i];
    //     let line = format!("{}\t{}\t{}\n", cur_taxon.0, cur_rinf, cur_err);
    //     file.write_all(line.as_bytes()).unwrap();
    // }

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