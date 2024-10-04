use crate::model::error::{DeltaPDError, DeltaPDResult};
use std::collections::HashMap;
use std::path::Path;

#[derive(Debug)]
pub struct MetadataFileRow {
    pub genome_id: String,
    pub contig_length: usize,
}

impl MetadataFileRow {
    pub fn new(genome_id: String, contig_length: usize) -> Self {
        MetadataFileRow {
            genome_id,
            contig_length,
        }
    }
}

pub struct MetadataFile {
    pub rows: HashMap<String, MetadataFileRow>,
}


impl MetadataFile {
    pub fn read(path: &Path, delimiter: u8) -> DeltaPDResult<Self> {
        // Create the output rows
        let mut rows: HashMap<String, MetadataFileRow> = HashMap::new();

        // Create the reader
        let reader = csv::ReaderBuilder::new()
            .delimiter(delimiter)
            .from_path(path)
            .map_err(DeltaPDError::CsvError)?;

        // Parse each line
        for (line_no, result) in reader.into_records().enumerate() {
            let record = result.map_err(DeltaPDError::CsvError)?;

            // Unique sequence identifier
            let sequence_id = record
                .get(0)
                .ok_or_else(|| {
                    DeltaPDError::Error(format!("No sequence id present at line {line_no}."))
                })?
                .to_string();

            // Genome identifier (may be non-unique)
            let genome_id = record
                .get(1)
                .ok_or_else(|| {
                    DeltaPDError::Error(format!("No genome id present at line {line_no}."))
                })?
                .to_string();

            // Contig length
            let contig_length = record
                .get(2)
                .ok_or_else(|| {
                    DeltaPDError::Error(format!(
                        "No contig length present at line {line_no}."
                    ))
                })?
                .parse()
                .map_err(|_| {
                    DeltaPDError::Error(format!(
                        "Could not parse contig length at line {line_no}."
                    ))
                })?;

            // Check for duplicate sequence IDs
            if rows.contains_key(&sequence_id) {
                return Err(DeltaPDError::Error(format!(
                    "Duplicate sequence ID found at line {line_no}."
                )));
            }

            // Store the row
            rows.insert(sequence_id, MetadataFileRow::new(genome_id, contig_length));
        }

        // Finally, return the metadata file
        Ok(MetadataFile {
            rows,
        })
    }

    pub fn get_ref_taxon(&self, qry_taxon: &str) -> Option<&str> {
        self.rows.get(qry_taxon).map(|row| row.genome_id.as_ref())
    }

    pub fn get_contig_len(&self, qry_taxon: &str) -> DeltaPDResult<usize> {
        self.rows.get(qry_taxon).map(|row| row.contig_length).ok_or_else(|| {
            DeltaPDError::Error(format!("No contig length found for taxon: {}", qry_taxon))
        })
    }
}
