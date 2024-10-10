use crate::model::linalg::{LinearModelCorr, LinearModelError, LinearModelType, PyLinearModelCorr, PyLinearModelError, PyLinearModelType};
use pyo3::{pyclass, pymethods};
use std::cmp;
use std::collections::HashSet;
use std::path::PathBuf;





pub struct Params {
    // Parameters
    pub cpus: usize,

    pub knn: usize,
    pub sample_size: f64,
    pub replicates: usize,
    pub taxa: HashSet<String>,

    // Model parameters
    pub model: LinearModelType,
    pub model_error: LinearModelError,
    pub model_corr: LinearModelCorr,
    pub direction: ParamsDirection,

    pub debug: bool,
    pub output_dir: PathBuf,
}


impl Params {
    pub fn new(
        cpus: usize,
        knn: usize,
        sample_size: f64,
        replicates: usize,
        taxa: HashSet<String>,
        model: LinearModelType,
        model_error: LinearModelError,
        model_corr: LinearModelCorr,
        direction: ParamsDirection,
        debug: bool,
        output_dir: PathBuf,
    ) -> Self {
        Self {
            cpus: cmp::max(cpus, 1),
            knn,
            sample_size,
            replicates: cmp::max(replicates, 1),
            taxa,
            model,
            model_error,
            model_corr,
            direction,
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
    pub fn new(
        cpus: usize,
        knn: usize,
        sample_size: f64,
        replicates: usize,
        taxa: HashSet<String>,
        model: PyLinearModelType,
        error: PyLinearModelError,
        corr: PyLinearModelCorr,
        direction: ParamsDirection,
        debug: bool,
        output_dir: String
    ) -> Self {
        PyParams {
            params: Params::new(
                cpus,
                knn,
                sample_size,
                replicates,
                taxa,
                model.to_enum(),
                error.to_enum(),
                corr.to_enum(),
                direction,
                debug,
                PathBuf::from(output_dir)
            )
        }
    }
}


#[derive(Clone)]
#[pyclass]
pub enum ParamsDirection {
    QvR,
    RvQ,
}
