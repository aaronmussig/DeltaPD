use crate::model::metadata::MetadataFile;
use crate::model::pdm::{QryDistMatrix, RefDistMatrix};
use bitvec::bitvec;
use bitvec::prelude::BitVec;
use phylodm::tree::Taxon;
use std::collections::{HashMap, HashSet};

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
    pub fn new(job: &'a Job, qry_mat: &'a QryDistMatrix, ref_mat: &'a RefDistMatrix, metadata_file: &MetadataFile, sample_size: f64) -> Self {

        // Extract the query taxa of interest
        let taxon = job.taxon;

        // Sample taxa with replacement
        let qry_taxa_sampled = qry_mat.dm.sample_taxa_with_replacement_include(sample_size, taxon);
        let qry_taxa_sampled_unq = qry_taxa_sampled.iter().cloned().collect::<HashSet<&Taxon>>();

        // Calculate the number of data points that will be included in this
        let n_data_points = qry_taxa_sampled.len() * (qry_taxa_sampled.len() - 1) / 2;

        // Collect the pairwise distances for those query taxa
        let mut qry_data: Vec<f64> = Vec::with_capacity(n_data_points);
        let mut ref_data: Vec<f64> = Vec::with_capacity(n_data_points);
        let mut qry_labels: Vec<(&Taxon, &Taxon)> = Vec::with_capacity(n_data_points);
        let mut ref_labels: Vec<(&Taxon, &Taxon)> = Vec::with_capacity(n_data_points);
        let mut label_bitmask: HashMap<&Taxon, BitVec> = HashMap::with_capacity(qry_taxa_sampled_unq.len());

        // Initialise the bitmask for each taxon
        for qry_taxon in qry_taxa_sampled_unq {
            label_bitmask.insert(qry_taxon, bitvec![0; n_data_points]);
        }

        // Iterate over each randomly sampled taxon index
        let mut cur_index = 0usize;
        for (idx_i, &qry_taxon_i) in qry_taxa_sampled.iter().enumerate() {

            // Collect the query and reference taxa (from)
            let ref_taxon_i = {
                let ref_taxon_i_str = metadata_file.get_ref_taxon(&qry_taxon_i.0).unwrap();
                ref_mat.dm.get_taxon_from_string(ref_taxon_i_str)
            };

            // Iterate over the lower triangle of the randomly sampled taxa indices
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
