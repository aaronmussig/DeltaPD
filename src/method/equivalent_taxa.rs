use crate::model::pdm::DistMatrix;
use phylodm::tree::Taxon;
use std::collections::HashMap;

pub fn find_equivalent_taxa(dist_matrix: &DistMatrix) -> HashMap<&Taxon, Vec<&Taxon>> {

    // This hasn't been tested but it should work!

    let mut taxon_to_equivalent_taxa: HashMap<&Taxon, Vec<&Taxon>> = HashMap::new();

    // Iterate over each taxon
    for (taxon_idx, taxon) in dist_matrix.taxa.iter().enumerate() {

        // Get the corresponding row in the distance matrix
        let row = dist_matrix.matrix.row(taxon_idx);

        // Find the index of any zero values
        let zero_idx = row.iter().position(|&x| x == 0.0);

        // Collect the taxa that correspond to those indices
        let equivalent_taxa: Vec<&Taxon> = zero_idx.iter().map(|&idx| &dist_matrix.taxa[idx]).collect();

        // Remove the current taxon from the list
        for eq_taxon in equivalent_taxa {

            // Ignore self hits
            if taxon == taxon {
                continue;
            }

            // Add the equivalent taxa to the hashmap
            taxon_to_equivalent_taxa.entry(taxon).or_insert(Vec::new()).push(eq_taxon);
        }
        //

    }
    taxon_to_equivalent_taxa
}