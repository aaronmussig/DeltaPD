use crate::model::error::{DeltaPDError, DeltaPDResult};
use csv::Writer;
use ndarray::Array2;
use phylodm::tree::{Edge, NodeId, Taxon};
use phylodm::PDM;
use pyo3::exceptions::PyValueError;
use pyo3::types::PyType;
use pyo3::{pyclass, pymethods, Bound, PyResult};
use rand::distributions::{Distribution, Uniform};
use rand::thread_rng;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use std::str::FromStr;
use std::time::Instant;

#[derive(Debug, Clone)]
pub struct DistMatrix {
    pub taxa: Vec<Taxon>,
    pub taxon_to_idx: HashMap<Taxon, usize>,
    pub taxon_str_to_idx: HashMap<String, usize>,
    pub matrix: Array2<f64>,
}

impl DistMatrix {
    pub fn new(taxa: Vec<Taxon>, matrix: Array2<f64>) -> DistMatrix {
        let mut taxon_to_idx: HashMap<Taxon, usize> = HashMap::with_capacity(taxa.len());
        let mut taxon_str_to_idx: HashMap<String, usize> = HashMap::with_capacity(taxa.len());

        for (i, taxon) in taxa.iter().enumerate() {
            taxon_to_idx.insert(taxon.clone(), i);
            taxon_str_to_idx.insert(taxon.0.clone(), i);
        }

        DistMatrix { taxa, taxon_to_idx, taxon_str_to_idx, matrix }
    }

    /// Samples sample_size number of taxa from the matrix with replacement. Always includes taxon.
    pub fn sample_taxa_with_replacement_include<'a>(&'a self, sample_size: f64, taxon: &'a Taxon) -> Vec<&'a Taxon> {
        let n = self.taxa.len();

        // If sample size is greater than 1, then we are taking an absolute number
        let n_samples = if sample_size > 1.0 {
            sample_size as usize
        } else {
            (n as f64 * sample_size) as usize
        };

        // Setup the random number generator
        let mut rng = thread_rng();
        let distribution = Uniform::new(0, n);

        // Create the output vector
        let mut out: Vec<&Taxon> = Vec::with_capacity(n_samples + 1);
        while out.len() < n_samples {
            let cur_taxon_idx = distribution.sample(&mut rng);
            let cur_taxon = &self.taxa[cur_taxon_idx];
            if cur_taxon != taxon {
                out.push(cur_taxon);
            }
        }
        // Add the taxon of interest
        out.push(taxon);
        out
    }

    pub fn get_taxon_from_string(&self, taxon: &str) -> &Taxon {
        let taxon_index = *self.taxon_str_to_idx.get(taxon).unwrap();
        &self.taxa[taxon_index]
    }

    pub fn distance(&self, a: &Taxon, b: &Taxon) -> f64 {
        let a_idx = *self.taxon_to_idx.get(a).unwrap();
        let b_idx = *self.taxon_to_idx.get(b).unwrap();
        self.matrix[[a_idx, b_idx]]
    }

    pub fn to_file(&self, path: &Path) -> DeltaPDResult<()> {
        let mut file = Writer::from_path(path).map_err(DeltaPDError::CsvError)?;

        // Write the header line (taxa names)
        file.write_field("").map_err(DeltaPDError::CsvError)?;
        for taxon in &self.taxa {
            file.write_field(taxon.0.as_str()).map_err(DeltaPDError::CsvError)?;
        }
        file.write_record(None::<&[u8]>).map_err(DeltaPDError::CsvError)?;

        // Write the matrix
        for (i, row) in self.matrix.outer_iter().enumerate() {
            file.write_field(self.taxa[i].0.as_str()).map_err(DeltaPDError::CsvError)?;
            for value in row.iter() {
                file.write_field(value.to_string()).map_err(DeltaPDError::CsvError)?;
            }
            file.write_record(None::<&[u8]>).map_err(DeltaPDError::CsvError)?;
        }
        Ok(())
    }

    pub fn from_path(path: &Path) -> DeltaPDResult<Self> {

        // Read the file
        let file = File::open(path).map_err(DeltaPDError::IoError)?;
        let reader = io::BufReader::new(file);

        let mut seen_taxa = false;
        let mut seen_edges = false;

        let mut file_taxa: HashMap<Taxon, usize> = HashMap::new();
        let mut file_edges: Vec<(usize, usize, f64)> = Vec::new();

        for line in reader.lines() {
            let line = line.map_err(DeltaPDError::IoError)?;

            // Check what section we are up to
            if line.starts_with("#TAXA") {
                seen_taxa = true;
                continue;
            } else if line.starts_with("#EDGES") {
                seen_edges = true;
                continue;
            } else {
                if seen_edges {
                    // We are adding the edges
                    let parts = line.split('\t').collect::<Vec<&str>>();
                    if parts.len() != 3 {
                        return Err(DeltaPDError::Error("Invalid number of columns in the edges section".to_string()));
                    }
                    let parent = usize::from_str(parts[0]).map_err(DeltaPDError::ParseIntError)?;
                    let child = usize::from_str(parts[1]).map_err(DeltaPDError::ParseIntError)?;
                    let length = f64::from_str(parts[2]).map_err(DeltaPDError::ParseFloatError)?;
                    file_edges.push((parent, child, length));
                } else if seen_taxa {
                    // We are adding the taxa
                    let parts = line.split('\t').collect::<Vec<&str>>();
                    if parts.len() != 2 {
                        return Err(DeltaPDError::Error("Invalid number of columns in the taxa section".to_string()));
                    }
                    let taxon = Taxon(parts[0].to_string());
                    let idx = usize::from_str(parts[1]).map_err(DeltaPDError::ParseIntError)?;
                    file_taxa.insert(taxon, idx);
                } else {
                    // This is an error as neither section has been found
                    return Err(DeltaPDError::Error("No section found in the file".to_string()));
                }
            }
        }

        // Create the PhyloDM object
        let mut pdm = PDM::default();

        // Record the old index to new index mapping
        let mut old_id_to_new_id: HashMap<usize, NodeId> = HashMap::with_capacity(file_edges.len());

        // Add the taxa
        for (taxon, idx) in file_taxa {
            let new_id = pdm.add_node(Some(&taxon)).unwrap();
            old_id_to_new_id.insert(idx, new_id);
        }

        // Add the edges
        for (parent_id, child_id, length) in file_edges {
            let parent_node_id = match old_id_to_new_id.get(&parent_id) {
                Some(id) => *id,
                None => {
                    let new_id = pdm.add_node(None).unwrap();
                    old_id_to_new_id.insert(parent_id, new_id);
                    new_id
                }
            };
            let child_node_id = match old_id_to_new_id.get(&child_id) {
                Some(id) => *id,
                None => {
                    let new_id = pdm.add_node(None).unwrap();
                    old_id_to_new_id.insert(child_id, new_id);
                    new_id
                }
            };
            pdm.add_edge(parent_node_id, child_node_id, Edge(length));
        }

        println!("Total sum of all branch lengths: {:?}", pdm.length());
        let previous_time = Instant::now();
        let (taxa, matrix) = pdm.matrix(true).unwrap();
        println!("Time to create PDM: {:?}", previous_time.elapsed());
        Ok(Self::new(taxa, matrix))
    }
}

#[pyclass]
#[derive(Debug, Clone)]
pub struct PyDistMatrix {
    pub dm: DistMatrix,
}

#[pymethods]
impl PyDistMatrix {
    #[new]
    pub fn new(taxa: Vec<(String, usize)>, edges: Vec<(usize, usize, f64)>) -> Self {
        let mut pdm = PDM::default();
        let mut old_id_to_new_id: HashMap<usize, NodeId> = HashMap::with_capacity(edges.len());

        // Create the named nodes
        for (taxon, id) in taxa {
            let taxon = Taxon(taxon);
            let new_id = pdm.add_node(Some(&taxon)).unwrap();
            old_id_to_new_id.insert(id, new_id);
        }

        // Create the unnamed nodes and edges
        for (parent_id, child_id, length) in edges {

            // Get or create the node id for the parent node
            let parent_node_id = match old_id_to_new_id.get(&parent_id) {
                Some(id) => *id,
                None => {
                    let new_id = pdm.add_node(None).unwrap();
                    old_id_to_new_id.insert(parent_id, new_id);
                    new_id
                }
            };

            // Get or create the node id for the child node
            let child_node_id = match old_id_to_new_id.get(&child_id) {
                Some(id) => *id,
                None => {
                    let new_id = pdm.add_node(None).unwrap();
                    old_id_to_new_id.insert(child_id, new_id);
                    new_id
                }
            };

            pdm.add_edge(parent_node_id, child_node_id, Edge(length));
        }

        println!("Total sum of all branch lengths: {:?}", pdm.length());

        let (dm_taxa, dm) = pdm.matrix(true).unwrap();
        let dm = DistMatrix::new(dm_taxa, dm);
        Self { dm }
    }

    pub fn to_file(&self, path: &str) -> PyResult<()> {
        self.dm.to_file(Path::new(path)).map_err(|e| PyValueError::new_err(format!("{:?}", e)))
    }

    #[classmethod]
    pub fn from_file(_cls: &Bound<'_, PyType>, path: &str) -> PyResult<Self> {
        let dm = DistMatrix::from_path(Path::new(path)).map_err(|e| PyValueError::new_err(format!("{:?}", e)))?;
        Ok(Self { dm })
    }

    #[classmethod]
    pub fn from_path(_cls: &Bound<'_, PyType>, path: &str) -> PyResult<Self> {
        let dm = DistMatrix::from_path(Path::new(path)).map_err(|e| PyValueError::new_err(format!("{:?}", e)))?;
        Ok(Self { dm })
    }
}


#[derive(Debug)]
pub struct QryDistMatrix {
    pub dm: DistMatrix,
}


#[derive(Debug)]
pub struct RefDistMatrix {
    pub dm: DistMatrix,
}

