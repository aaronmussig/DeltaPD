use crate::method::knn::{KnnVecAsLinalg, PyKnnVecAsLinalg};
use crate::model::linalg::{LinearModel, PyLinearModel};
use phylodm::tree::Taxon;
use pyo3::pyclass;

pub struct OutputResult {
    pub linear_model_base: LinearModel,
    pub linear_model_jk: Vec<LinearModel>,
    pub relative_influence: Vec<f64>,
    pub std_error: Vec<f64>,
    pub query_taxon: Taxon,
    pub knn_as_linalg: KnnVecAsLinalg,
}


impl OutputResult {
    pub fn to_python(&self) -> PyOutputResult {
        PyOutputResult {
            linear_model_base: self.linear_model_base.to_python(),
            linear_model_jk: self.linear_model_jk.iter().map(|x| x.to_python()).collect(),
            relative_influence: self.relative_influence.clone(),
            std_error: self.std_error.clone(),
            query_taxon: self.query_taxon.0.clone(),
            knn_as_linalg: self.knn_as_linalg.to_python(),
        }
    }
}


#[pyclass]
pub struct PyOutputResult {
    #[pyo3(get, set)]
    pub linear_model_base: PyLinearModel,

    #[pyo3(get, set)]
    pub linear_model_jk: Vec<PyLinearModel>,

    #[pyo3(get, set)]
    pub relative_influence: Vec<f64>,

    #[pyo3(get, set)]
    pub std_error: Vec<f64>,

    #[pyo3(get, set)]
    pub query_taxon: String,

    #[pyo3(get, set)]
    pub knn_as_linalg: PyKnnVecAsLinalg,
}


pub struct OutputResultSmall {
    pub taxa: Vec<Taxon>,
    pub std_error: Vec<f64>,
}


impl OutputResultSmall {
    pub fn to_python(&self) -> PyOutputResultSmall {
        PyOutputResultSmall {
            std_error: self.std_error.clone(),
            taxa: self.taxa.iter().map(|x| x.0.clone()).collect(),
        }
    }
}


#[pyclass]
pub struct PyOutputResultSmall {
    #[pyo3(get, set)]
    pub std_error: Vec<f64>,

    #[pyo3(get, set)]
    pub taxa: Vec<String>,

}