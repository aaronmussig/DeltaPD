use clap::ValueEnum;
use pyo3::pyclass;

use crate::model::types::{CorrFn, ErrorFn, ModelFn};
use crate::stats::linalg::{
    calc_mse, calc_mse_norm, calc_pearson_correlation, calc_r2, calc_repeated_median, calc_rmse,
    calc_theil_sen_gradient,
};

#[derive(Copy, Clone, Debug, ValueEnum)]
pub enum LinearModelType {
    /// Repeated median
    RepeatedMedian,
    /// Theil-Sen
    TheilSen,
}


impl LinearModelType {
    pub fn get_fn(&self) -> ModelFn {
        match self {
            LinearModelType::RepeatedMedian => calc_repeated_median,
            LinearModelType::TheilSen => calc_theil_sen_gradient,
        }
    }
}

#[derive(Clone)]
#[pyclass]
pub enum PyLinearModelType {
    RepeatedMedian,
    TheilSen,
}

impl PyLinearModelType {
    pub fn to_enum(&self) -> LinearModelType {
        match self {
            PyLinearModelType::RepeatedMedian => LinearModelType::RepeatedMedian,
            PyLinearModelType::TheilSen => LinearModelType::TheilSen,
        }
    }
}


#[derive(Copy, Clone, ValueEnum)]
pub enum LinearModelError {
    MSE,
    NormMSE,
    RMSE,
}


impl LinearModelError {
    pub fn get_fn(&self) -> ErrorFn {
        match self {
            LinearModelError::MSE => calc_mse,
            LinearModelError::NormMSE => calc_mse_norm,
            LinearModelError::RMSE => calc_rmse,
        }
    }
}

#[derive(Clone)]
#[pyclass]
pub enum PyLinearModelError {
    MSE,
    NormMSE,
    RMSE,
}

impl PyLinearModelError {
    pub fn to_enum(&self) -> LinearModelError {
        match self {
            PyLinearModelError::MSE => LinearModelError::MSE,
            PyLinearModelError::NormMSE => LinearModelError::NormMSE,
            PyLinearModelError::RMSE => LinearModelError::RMSE,
        }
    }
}


#[derive(Copy, Clone, ValueEnum)]
pub enum LinearModelCorr {
    R2,
    Pearson,
}


impl LinearModelCorr {
    pub fn get_fn(&self) -> CorrFn {
        match self {
            LinearModelCorr::R2 => calc_r2,
            LinearModelCorr::Pearson => calc_pearson_correlation,
        }
    }
}


#[derive(Clone)]
#[pyclass]
pub enum PyLinearModelCorr {
    R2,
    Pearson,
}

impl PyLinearModelCorr {
    pub fn to_enum(&self) -> LinearModelCorr {
        match self {
            PyLinearModelCorr::R2 => LinearModelCorr::R2,
            PyLinearModelCorr::Pearson => LinearModelCorr::Pearson,
        }
    }
}

pub struct LinearModelParams {
    pub gradient: f64,
    pub intercept: f64,
}

impl LinearModelParams {
    pub fn new(gradient: f64, intercept: f64) -> Self {
        Self {
            gradient,
            intercept,
        }
    }

    pub fn to_python(&self) -> PyLinearModelParams {
        PyLinearModelParams {
            gradient: self.gradient,
            intercept: self.intercept
        }
    }
}

#[derive(Clone)]
#[pyclass]
pub struct PyLinearModelParams {
    #[pyo3(get, set)]
    pub gradient: f64,
    #[pyo3(get, set)]
    pub intercept: f64
}

pub struct LinearModelEval {
    pub error: f64,
    pub corr: f64,
}

impl LinearModelEval {
    pub fn new(error: f64, corr: f64) -> Self {
        Self { error, corr }
    }

    pub fn to_python(&self) -> PyLinearModelEval {
        PyLinearModelEval {
            error: self.error,
            corr: self.corr
        }
    }
}

#[derive(Clone)]
#[pyclass]
pub struct PyLinearModelEval {
    #[pyo3(get, set)]
    pub error: f64,
    #[pyo3(get, set)]
    pub corr: f64
}


pub struct LinearModel {
    pub params: LinearModelParams,
    pub eval: LinearModelEval,
}

impl LinearModel {
    pub fn fit(
        model_type: LinearModelType,
        model_error: LinearModelError,
        model_corr: LinearModelCorr,
        x: &[f64],
        y: &[f64],
    ) -> Self {
        // Fit a linear model
        let model_fn = model_type.get_fn();
        let params = model_fn(x, y);

        // Calculate the error
        let error_fn = model_error.get_fn();
        let error = error_fn(x, y);

        // Calculate the correlation
        let corr_fn = model_corr.get_fn();
        let corr = corr_fn(x, y);

        Self {
            params,
            eval: LinearModelEval { error, corr },
        }
    }

    pub fn predict(&self, x: &[f64]) -> Vec<f64> {
        x.iter()
            .map(|x0| self.params.gradient * x0 + self.params.intercept)
            .collect()
    }

    pub fn to_python(&self) -> PyLinearModel {
        PyLinearModel {
            params: self.params.to_python(),
            eval: self.eval.to_python()
        }
    }
}

#[derive(Clone)]
#[pyclass]
pub struct PyLinearModel {
    #[pyo3(get, set)]
    pub params: PyLinearModelParams,
    #[pyo3(get, set)]
    pub eval: PyLinearModelEval
}

