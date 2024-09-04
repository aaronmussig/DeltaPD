use std::cmp;
use std::collections::{HashMap, HashSet};

use log::warn;
use ndarray::Array2;
use phylodm::tree::Taxon;
use pyo3::pyclass;
use crate::model::error::DeltaPDResult;
use crate::model::metadata::MetadataFile;
use crate::model::pdm::{QryDistMatrix, RefDistMatrix};
use crate::ndarray::sort::get_nn_from_distance_matrix;

use rand::thread_rng;
use rand::seq::SliceRandom;

/// Get the nearest neighbours of all leaf nodes by patristic distance in a PDM.
pub fn get_pdm_k_bootstrap_from_matrix(
    taxa: &[Taxon],
    limit: &HashSet<&Taxon>,
    matrix: &Array2<f64>,
    k: usize,
) -> DeltaPDResult<HashMap<Taxon, Vec<usize>>> {
    // Sanity check
    if k > taxa.len() {
        warn!("The value of k-nearest neighbours is greater than the number of taxa in the PDM. The full set of taxa will be used.");
    }
    let knn = cmp::min(k, taxa.len());

    println!("Getting bootstraps");
    // Maybe add a distribution for sampling based on the nearest neighbours (like a skew right).


    for (taxon_idx, taxon) in taxa.iter().enumerate() {

        // Get the corresponding slice of the matrix for this taxon
        let row = matrix.row(taxon_idx);




    }


    // Obtain the indices for the nearest neighbours for all taxa
    let taxon_idx_to_nn = get_nn_from_distance_matrix(&matrix.view());

    // Generate the output data structure
    let mut out: HashMap<Taxon, Vec<usize>> = HashMap::with_capacity(taxa.len());

    // Iterate over each taxon and the nearest neighbours
    for (from_idx, taxon_vec) in taxon_idx_to_nn.iter().enumerate() {
        // This is the taxon we are looking at, we will now iterate over its nearest neighbours
        let cur_outer_taxon = &taxa[from_idx];

        // Skip processing this taxon if it's not in the list
        if !limit.contains(cur_outer_taxon) {
            continue;
        }

        // Generate the list of nearest neighbours
        let mut cur_row: Vec<usize> = Vec::with_capacity(knn);

        // Shuffle taxon_vec
        let mut taxon_vec = taxon_vec.clone();
        taxon_vec.shuffle(&mut thread_rng());

        for taxon_idx in &taxon_vec {
            // Check that this taxon is in the allowed set
            let cur_inner_taxon = &taxa[*taxon_idx];
            if !limit.contains(cur_inner_taxon) {
                continue;
            }

            // Stop once we have obtained k taxa
            if cur_row.len() >= knn {
                break;
            }

            // Don't keep self-comparisons (always 0)
            if *taxon_idx != from_idx {
                cur_row.push(*taxon_idx);
            }
        }
        out.insert(cur_outer_taxon.clone(), cur_row);
    }
    Ok(out)
}



/// Get the nearest neighbours of all leaf nodes by patristic distance in a PDM.
pub fn get_pdm_k_nearest_neighbours_from_matrix(
    taxa: &[Taxon],
    limit: &HashSet<&Taxon>,
    matrix: &Array2<f64>,
    k: usize,
) -> DeltaPDResult<HashMap<Taxon, Vec<usize>>> {
    // Sanity check
    if k > taxa.len() {
        warn!("The value of k-nearest neighbours is greater than the number of taxa in the PDM. The full set of taxa will be used.");
    }
    let knn = cmp::min(k, taxa.len());

    // Obtain the indices for the nearest neighbours for all taxa
    let taxon_idx_to_nn = get_nn_from_distance_matrix(&matrix.view());

    // Generate the output data structure
    let mut out: HashMap<Taxon, Vec<usize>> = HashMap::with_capacity(taxa.len());

    // Iterate over each taxon and the nearest neighbours
    for (from_idx, taxon_vec) in taxon_idx_to_nn.iter().enumerate() {
        // This is the taxon we are looking at, we will now iterate over its nearest neighbours
        let cur_outer_taxon = &taxa[from_idx];

        // Skip processing this taxon if it's not in the list
        if !limit.contains(cur_outer_taxon) {
            continue;
        }

        // Generate the list of nearest neighbours
        let mut cur_row: Vec<usize> = Vec::with_capacity(knn);
        for taxon_idx in taxon_vec {
            // Check that this taxon is in the allowed set
            let cur_inner_taxon = &taxa[*taxon_idx];
            if !limit.contains(cur_inner_taxon) {
                continue;
            }

            // Stop once we have obtained k taxa
            if cur_row.len() >= knn {
                break;
            }

            // Don't keep self-comparisons (always 0)
            if *taxon_idx != from_idx {
                cur_row.push(*taxon_idx);
            }
        }
        out.insert(cur_outer_taxon.clone(), cur_row);
    }
    Ok(out)
}


pub struct KnnVecAsLinalg {
    pub qry_data: Vec<f64>,
    pub ref_data: Vec<f64>,
    pub qry_labels: Vec<Taxon>,
    pub ref_labels: Vec<Taxon>
}

impl KnnVecAsLinalg {
    pub fn to_python(&self) -> PyKnnVecAsLinalg {
        PyKnnVecAsLinalg {
            qry_data: self.qry_data.clone(),
            ref_data: self.ref_data.clone(),
            qry_labels: self.qry_labels.iter().map(|x| x.0.clone()).collect(),
            ref_labels: self.ref_labels.iter().map(|x| x.0.clone()).collect()
        }
    }
}

#[derive(Clone)]
#[pyclass]
pub struct PyKnnVecAsLinalg {
    #[pyo3(get, set)]
    pub qry_data: Vec<f64>,
    #[pyo3(get, set)]
    pub ref_data: Vec<f64>,
    #[pyo3(get, set)]
    pub qry_labels: Vec<String>,
    #[pyo3(get, set)]
    pub ref_labels: Vec<String>

}



/// For a given taxon
pub fn create_vecs_from_knn2(
    qry_mat: &QryDistMatrix,
    ref_mat: &RefDistMatrix,
    q_taxon: &Taxon,
    metadata_file: &MetadataFile,
    knn_qry: &HashMap<Taxon, Vec<usize>>,
) -> DeltaPDResult<KnnVecAsLinalg> {
    // For the given query taxon, obtain the KNN vector of indices
    let q_knn_vec = knn_qry.get(q_taxon).unwrap();
    let n_knn = q_knn_vec.len();

    // Create the output vectors
    let mut qry_vec: Vec<f64> = Vec::with_capacity(n_knn);
    let mut ref_vec: Vec<f64> = Vec::with_capacity(n_knn);
    let mut qry_labels: Vec<Taxon> = Vec::with_capacity(n_knn);
    let mut ref_labels: Vec<Taxon> = Vec::with_capacity(n_knn);

    // Get the reference taxon since this is what we will be comparing to
    let r_taxon = metadata_file.get_ref_taxon(&q_taxon.0).unwrap();

    // For the given query taxon, obtain the KNN vector of indices
    let q_knn_vec = knn_qry.get(q_taxon).unwrap();

    // Iterate over each index to obtain the corresponding taxon
    for q_knn_idx in q_knn_vec {
        // Obtain the taxon from the index
        let q_knn_taxon = &qry_mat.dm.taxa[*q_knn_idx];

        // Obtain the corresponding taxon from the reference PDM
        let r_knn_taxon = metadata_file.get_ref_taxon(&q_knn_taxon.0).unwrap();

        // Obtain the distance between the query taxon and the KNN taxon
        let q_dist = qry_mat.dm.distance(&q_taxon, &q_knn_taxon);

        // Obtain the distance between the reference taxon and the KNN taxon
        let r_dist = ref_mat.dm.distance(&Taxon(r_taxon.to_string()), &Taxon(r_knn_taxon.to_string()));

        // Push the distances to the vectors
        qry_vec.push(q_dist);
        ref_vec.push(r_dist);

        // Push the labels to the vector
        qry_labels.push(q_knn_taxon.clone());
        ref_labels.push(Taxon(r_knn_taxon.to_string()));
    }

    Ok(KnnVecAsLinalg {
        qry_data: qry_vec,
        ref_data: ref_vec,
        qry_labels,
        ref_labels
    })
}





