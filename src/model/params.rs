use pyo3::{pyclass, pymethods};
use crate::model::linalg::{LinearModelCorr, LinearModelError, LinearModelType, PyLinearModelCorr, PyLinearModelError, PyLinearModelType};

pub struct Params {
    // Parameters
    pub knn: usize,
    pub cpus: usize,

    // Model parameters
    pub model: LinearModelType,
    pub model_error: LinearModelError,
    pub model_corr: LinearModelCorr,
}


impl Params {
    pub fn new(
        knn: usize,
        cpus: usize,
        model: LinearModelType,
        model_error: LinearModelError,
        model_corr: LinearModelCorr,
    ) -> Self {
        Self {
            knn,
            cpus,
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
    pub fn new(knn: usize, cpus: usize, model: PyLinearModelType, error: PyLinearModelError, corr: PyLinearModelCorr) -> Self {
        PyParams {
            params: Params::new(knn, cpus, model.to_enum(), error.to_enum(), corr.to_enum())
        }
    }

}

