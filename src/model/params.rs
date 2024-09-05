use pyo3::{pyclass, pymethods};
use crate::model::linalg::{LinearModelCorr, LinearModelError, LinearModelType, PyLinearModelCorr, PyLinearModelError, PyLinearModelType};

pub struct Params {
    // Parameters
    pub cpus: usize,

    pub sample_size: f64,
    pub replicates: usize,
    pub taxa:Vec<String>,

    // Model parameters
    pub model: LinearModelType,
    pub model_error: LinearModelError,
    pub model_corr: LinearModelCorr,
}


impl Params {
    pub fn new(
        cpus: usize,
        sample_size: f64,
        replicates: usize,
        taxa: Vec<String>,
        model: LinearModelType,
        model_error: LinearModelError,
        model_corr: LinearModelCorr,
    ) -> Self {
        Self {
            cpus,
            sample_size,
            replicates,
            taxa,
            model,
            model_error,
            model_corr,
        }
    }
}

#[pyclass]
pub struct PyParams {
    pub params: Params
}

#[pymethods]
impl PyParams {

    #[new]
    pub fn new(cpus: usize, sample_size: f64, replicates: usize, taxa: Vec<String>, model: PyLinearModelType, error: PyLinearModelError, corr: PyLinearModelCorr) -> Self {
        PyParams {
            params: Params::new(cpus, sample_size, replicates, taxa, model.to_enum(), error.to_enum(), corr.to_enum())
        }
    }

}

