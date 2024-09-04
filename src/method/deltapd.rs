use std::collections::HashMap;
use std::time::Instant;
use phylodm::tree::Taxon;
use crate::method::common_taxa::find_common_ids_in_pdms;
use crate::method::knn::{create_vecs_from_knn2, get_pdm_k_bootstrap_from_matrix, get_pdm_k_nearest_neighbours_from_matrix};
use crate::model::error::DeltaPDResult;
use crate::model::linalg::LinearModel;
use crate::model::metadata::MetadataFile;
use crate::model::params::Params;
use crate::model::pdm::{QryDistMatrix, RefDistMatrix};
use crate::model::result::OutputResult;
use crate::stats::jackknife::{calc_prop_std_error_of_jackknife, calc_relative_jackknife_influence, jackknife_fn};
use rayon::prelude::*;

pub fn run_deltapd(qry_mat: &QryDistMatrix, ref_mat: &RefDistMatrix, metadata_file: &MetadataFile, params: &Params) -> DeltaPDResult<Vec<OutputResult>> {

    // Identify the common taxa between both trees
    let previous_time = Instant::now();
    let common_ids = find_common_ids_in_pdms(&ref_mat.dm.taxa, &qry_mat.dm.taxa, &metadata_file).unwrap();
    println!("Found common IDs in {:?}", previous_time.elapsed());

    // So far no cases have been found where the equivalent taxa are needed
    // let previous_time = Instant::now();
    // let equivalent_ref_taxa = find_equivalent_taxa(&ref_mat.dm);
    // let equivalent_qry_taxa = find_equivalent_taxa(&qry_mat.dm);
    // println!("Found equivalent taxa in {:?}", previous_time.elapsed());

    // Get the knn
    let previous_time = Instant::now();
    let knn_qry = get_pdm_k_bootstrap_from_matrix(
        &qry_mat.dm.taxa,
        &common_ids.qry_taxa,
        &qry_mat.dm.matrix,
        params.knn,
    ).unwrap();
    println!("Found KNNs in {:?}", previous_time.elapsed());

    println!("Iterating over each taxon (par)");
    let previous_time = Instant::now();
    let output: Vec<OutputResult> = knn_qry.par_iter().map(|(q_taxon, _qry_knn_idx)| {
        run_deltapd_on_taxon(q_taxon, qry_mat, ref_mat, metadata_file, params, &knn_qry).unwrap()
    }).collect();
    println!("Ran main method (par) in {:?}", previous_time.elapsed());

    Ok(output)
}

pub fn run_deltapd_on_taxon(q_taxon: &Taxon, qry_mat: &QryDistMatrix, ref_mat: &RefDistMatrix, metadata_file: &MetadataFile, params: &Params, knn_qry: &HashMap<Taxon, Vec<usize>>) -> DeltaPDResult<OutputResult> {

    // Create the XY vectors containing the data points from the KNN set
    let knn_as_linalg = create_vecs_from_knn2(&qry_mat, &ref_mat, q_taxon, metadata_file, &knn_qry)?;

    // Fit a model without any jackknifing
    let linear_model_base = LinearModel::fit(
        params.model,
        params.model_error,
        params.model_corr,
        &knn_as_linalg.ref_data,
        &knn_as_linalg.qry_data,
    );

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

    Ok(OutputResult {
        linear_model_base,
        linear_model_jk: jackknife_models,
        relative_influence,
        std_error,
        query_taxon: q_taxon.clone(),
        knn_as_linalg,
    })

}