use std::collections::HashMap;
use std::fs::File;
use std::path::Path;
use std::str::FromStr;

use csv::Writer;
use ndarray::Array2;
use phylodm::PDM;
use phylodm::tree::{Edge, NodeId, Taxon};
use pyo3::{Bound, pyclass, pymethods, PyResult};
use pyo3::exceptions::PyValueError;
use pyo3::types::PyType;

use crate::model::error::{DeltaPDError, DeltaPDResult};

#[derive(Debug, Clone)]
pub struct DistMatrix {
    pub taxa: Vec<Taxon>,
    pub taxon_to_idx: HashMap<Taxon, usize>,
    pub matrix: Array2<f64>,
}

impl DistMatrix {
    pub fn new(taxa: Vec<Taxon>, matrix: Array2<f64>) -> DistMatrix {
        let taxon_to_idx = taxa.iter().enumerate().map(|(i, taxon)| (taxon.clone(), i)).collect();
        DistMatrix { taxa, matrix, taxon_to_idx }
    }

    pub fn distance(&self, a: &Taxon, b: &Taxon) -> f64 {
        let a_idx = self.taxon_to_idx.get(a).unwrap();
        let b_idx = self.taxon_to_idx.get(b).unwrap();
        self.matrix[[*a_idx, *b_idx]]
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
        let mut reader = csv::Reader::from_reader(file);

        // Parse the header line
        let taxa: Vec<Taxon> = reader.headers().map_err(DeltaPDError::CsvError)?
            .iter()
            .skip(1)
            .map(|name| Taxon(name.to_string()))
            .collect();

        // Create the distance matrix
        let mut matrix = Array2::zeros((taxa.len(), taxa.len()));
        for (i, result) in reader.records().enumerate() {
            let records = result.map_err(DeltaPDError::CsvError)?;
            for (j, record) in records.iter().enumerate() {
                // The first j is the taxon name again
                if j > 0 {
                    let value = f64::from_str(record).map_err(DeltaPDError::ParseFloatError)?;
                    matrix[[i, j - 1]] = value;
                    matrix[[j - 1, i]] = value;
                }
            }
        }

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

