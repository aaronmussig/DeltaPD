use crate::model::metadata::MetadataFile;
use crate::model::pdm::{QryDistMatrix, RefDistMatrix};
use bitvec::bitvec;
use bitvec::prelude::BitVec;
use phylodm::tree::Taxon;
use std::collections::{HashMap, HashSet};
use crate::model::params::{Params, ParamsDirection};

pub struct Job<'a> {
    pub taxon: &'a Taxon,
    pub replicate: usize,
}

impl<'a> Job<'a> {
    pub fn new(taxon: &'a Taxon, replicate: usize) -> Self {
        Job {
            taxon,
            replicate,
        }
    }
}


pub struct JobData<'a> {
    pub taxon: &'a Taxon,
    pub ref_data: Vec<f64>,
    pub qry_data: Vec<f64>,
    pub qry_labels: Vec<(&'a Taxon, &'a Taxon)>,
    pub ref_labels: Vec<(&'a Taxon, &'a Taxon)>,
    pub label_bitmask: HashMap<&'a Taxon, BitVec>,
}

impl<'a> JobData<'a> {

    pub fn new(job: &'a Job, qry_mat: &'a QryDistMatrix, ref_mat: &'a RefDistMatrix, metadata_file: &MetadataFile, params: &Params) -> Self {

        // Extract the query taxa of interest
        let taxon = job.taxon;
        let taxon_ref_str = metadata_file.get_ref_taxon(&taxon.0).unwrap();
        let taxon_ref_obj = ref_mat.dm.get_taxon_from_string(taxon_ref_str);

        // Get the query taxa for the k-nearest neighbours
        let qry_taxa_sampled = match params.direction {

            // Query vs reference
            ParamsDirection::QvR => {

                // Get the closest query taxa to the query taxon
                let nearest_taxa = qry_mat.dm.argsort_taxon(taxon);

                // Keep track of what reference taxa we have seen, as we don't want to include
                // any query taxa that share the same reference taxon
                let mut ref_taxa_seen: HashSet<&str> = HashSet::new();

                // Get the k-nearest neighbours to the query taxon that do not share the same reference taxon
                let mut target_query_taxa: Vec<&Taxon> = Vec::with_capacity(params.knn);

                // Insert the current query taxon
                target_query_taxa.push(taxon);
                ref_taxa_seen.insert(taxon_ref_str);

                // Go over each of the nearest taxa
                for nearest_taxon in nearest_taxa {

                    // Find the reference taxon for the given query taxon
                    let cur_taxon_ref_str = metadata_file.get_ref_taxon(&nearest_taxon.0).unwrap();

                    // If we haven't seen the reference taxon yet include it
                    if !ref_taxa_seen.contains(cur_taxon_ref_str) {
                        target_query_taxa.push(nearest_taxon);
                        ref_taxa_seen.insert(cur_taxon_ref_str);
                    }

                    // Stop collecting query taxa if we've reached the target number of neighbours
                    if target_query_taxa.len() >= params.knn {
                        break;
                    }
                }
                target_query_taxa
            }

            // Reference vs query
            ParamsDirection::RvQ => {

                // TODO: In this case it would probably be good to exclude throwing in too many
                // TODO: query taxa that belong to the same reference taxon
                // TODO: Providing that the query taxa are sufficiently close

                // TODO: Take the nearest query taxa that belong to the reference taxa?

                let nearest_taxa = ref_mat.dm.argsort_taxon(taxon_ref_obj);

                // Get the k-nearest neighbours to the reference taxon
                let mut target_query_taxa: Vec<&Taxon> = Vec::with_capacity(params.knn);
                let mut n_ref_taxa_added = 0;

                for nearest_taxon in nearest_taxa {
                    let ref_qry_taxa = metadata_file.get_qry_taxa(&nearest_taxon.0);
                    if let Some(ref_qry_taxa) = ref_qry_taxa {
                        for qry_taxon_str in ref_qry_taxa {
                            let qry_taxon = qry_mat.dm.get_taxon_from_string(qry_taxon_str);
                            target_query_taxa.push(qry_taxon);
                        }
                        n_ref_taxa_added += 1;
                    }
                    if n_ref_taxa_added >= params.knn {
                        break;
                    }
                }
                target_query_taxa
            }
        };

        // Depending on the direction, get the k-nearest neighbors to this query taxon
        // qvr = getting the k nearest neighbours to the query taxon, omitting ones that share the same reference taxon
        // rvq = get the k nearest neighbours to the reference taxon then throwing in the query taxa that belong to it

        // Sample taxa with replacement
        // let qry_taxa_sampled = qry_mat.dm.sample_taxa_with_replacement_include(params.sample_size, taxon);
        let qry_taxa_sampled_unq = qry_taxa_sampled.iter().cloned().collect::<HashSet<&Taxon>>();

        // Calculate the number of data points that will be included in this
        let n_data_points = qry_taxa_sampled.len() * (qry_taxa_sampled.len() - 1) / 2;

        // Collect the pairwise distances for those query taxa
        let mut qry_data: Vec<f64> = Vec::with_capacity(n_data_points);
        let mut ref_data: Vec<f64> = Vec::with_capacity(n_data_points);
        let mut qry_labels: Vec<(&Taxon, &Taxon)> = Vec::with_capacity(n_data_points);
        let mut ref_labels: Vec<(&Taxon, &Taxon)> = Vec::with_capacity(n_data_points);
        let mut label_bitmask: HashMap<&Taxon, BitVec> = HashMap::with_capacity(qry_taxa_sampled_unq.len());

        // Initialise the bitmask for each taxon, i.e. which data points are associated with it
        for qry_taxon in qry_taxa_sampled_unq {
            label_bitmask.insert(qry_taxon, bitvec![0; n_data_points]);
        }

        // Iterate over each query taxon index
        let mut cur_index = 0usize;
        for (idx_i, &qry_taxon_i) in qry_taxa_sampled.iter().enumerate() {

            // Get the reference taxon associated with this query taxon (from)
            let ref_taxon_i = {
                let ref_taxon_i_str = metadata_file.get_ref_taxon(&qry_taxon_i.0).unwrap();
                ref_mat.dm.get_taxon_from_string(ref_taxon_i_str)
            };

            // Iterate over the lower triangle of the randomly sampled taxa indices matrix
            for idx_j in 0..idx_i {

                // Collect the query and reference taxa (to)
                let qry_taxon_j = qry_taxa_sampled[idx_j];
                let ref_taxon_j = {
                    let ref_taxon_j_str = metadata_file.get_ref_taxon(&qry_taxon_j.0).unwrap();
                    ref_mat.dm.get_taxon_from_string(ref_taxon_j_str)
                };

                // Get the distance between the two taxa
                let qry_dist = qry_mat.dm.distance(qry_taxon_i, qry_taxon_j);
                let ref_dist = ref_mat.dm.distance(ref_taxon_i, ref_taxon_j);

                // Store the values
                qry_data.push(qry_dist);
                ref_data.push(ref_dist);
                qry_labels.push((qry_taxon_i, qry_taxon_j));
                ref_labels.push((ref_taxon_i, ref_taxon_j));

                // Set the bitmask
                label_bitmask.get_mut(qry_taxon_i).unwrap().set(cur_index, true);
                label_bitmask.get_mut(qry_taxon_j).unwrap().set(cur_index, true);

                cur_index += 1;
            }
        }

        JobData {
            taxon,
            ref_data,
            qry_data,
            qry_labels,
            ref_labels,
            label_bitmask,
        }
    }
}
