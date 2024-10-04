use crate::model::linalg::{LinearModelCorr, LinearModelError, LinearModelType, PyLinearModelCorr, PyLinearModelError, PyLinearModelType};
use pyo3::{pyclass, pymethods};
use std::cmp;
use std::collections::HashSet;
use std::path::PathBuf;

pub struct Params {
    // Parameters
    pub cpus: usize,

    pub sample_size: f64,
    pub replicates: usize,
    pub taxa: HashSet<String>,

    // Model parameters
    pub model: LinearModelType,
    pub model_error: LinearModelError,
    pub model_corr: LinearModelCorr,

    pub debug: bool,
    pub output_dir: PathBuf,
}


impl Params {
    pub fn new(
        cpus: usize,
        sample_size: f64,
        replicates: usize,
        taxa: HashSet<String>,
        model: LinearModelType,
        model_error: LinearModelError,
        model_corr: LinearModelCorr,
        debug: bool,
        output_dir: PathBuf,
    ) -> Self {
        Self {
            cpus: cmp::max(cpus, 1),
            sample_size,
            replicates: cmp::max(replicates, 1),
            taxa,
            model,
            model_error,
            model_corr,
            debug,
            output_dir,
        }
    }
}

#[pyclass]
pub struct PyParams {
    pub params: Params,
}

#[pymethods]
impl PyParams {
    #[new]
    pub fn new(cpus: usize, sample_size: f64, replicates: usize, taxa: HashSet<String>, model: PyLinearModelType, error: PyLinearModelError, corr: PyLinearModelCorr, debug: bool, output_dir: String) -> Self {
        PyParams {
            params: Params::new(cpus, sample_size, replicates, taxa, model.to_enum(), error.to_enum(), corr.to_enum(), debug, PathBuf::from(output_dir))
        }
    }
}

