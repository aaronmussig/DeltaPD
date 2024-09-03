use std::collections::{HashMap, HashSet};

use log::{error, warn};
use phylodm::tree::Taxon;

use crate::model::error::DeltaPDResult;
use crate::model::metadata::MetadataFile;

pub struct CommonTaxa<'a> {
    pub ref_taxa: HashSet<&'a Taxon>,
    pub qry_taxa: HashSet<&'a Taxon>,
}

/// This method takes two vectors of taxa, and uses the mapping provided in the metadata
/// file to find those common.
pub fn find_common_ids_in_pdms<'a>(
    r_taxa: &'a [Taxon],
    q_taxa: &'a [Taxon],
    metadata_file: &MetadataFile,
) -> DeltaPDResult<CommonTaxa<'a>> {

    // Create a hashmap to get the pointer to the reference taxon
    let ref_taxa_map = {
        let mut out = HashMap::new();
        for taxon in r_taxa {
            out.insert(taxon.0.clone(), taxon);
        }
        out
    };

    // Create the output sets
    let mut out_ref: HashSet<&'a Taxon> = HashSet::with_capacity(r_taxa.len());
    let mut out_qry: HashSet<&'a Taxon> = HashSet::with_capacity(q_taxa.len());

    // Iterate over the query taxa and match their corresponding reference taxon
    for qry_taxon in q_taxa {
        // Read the reference ID from the metadata file
        if let Some(ref_taxon_str) = metadata_file.get_ref_taxon(&qry_taxon.0) {
            // Match the reference taxon to that present in the reference taxa set
            if let Some(ref_taxon) = ref_taxa_map.get(ref_taxon_str) {
                out_ref.insert(ref_taxon);
                out_qry.insert(qry_taxon);
            } else {
                error!(
                    "Reference taxon '{}' not found in the reference tree, skipping.",
                    ref_taxon_str
                );
            }
        } else {
            warn!(
                "Query taxon '{}' has no reference taxon set in the metadata file, skipping.",
                qry_taxon.0
            );
        }
    }

    // Return
    Ok(CommonTaxa {
        ref_taxa: out_ref,
        qry_taxa: out_qry,
    })
}

