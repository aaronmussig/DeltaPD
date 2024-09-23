use crate::model::error::{DeltaPDError, DeltaPDResult};
use crate::model::types::{CorrFn, ErrorFn, ModelFn};
use crate::ndarray::sort::argsort_by_vec;
use crate::stats::linalg::{calc_mse, calc_mse_norm, calc_pearson_correlation, calc_r2, calc_repeated_median, calc_rmse, calc_theil_sen_gradient};
use crate::stats::vec::{calc_median, calc_median_sorted};
use crate::util::bitvec::bitvec_boolean_a_and_not_b;
use bitvec::bitvec;
use bitvec::vec::BitVec;
use clap::ValueEnum;
use pyo3::pyclass;

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

#[derive(Clone, Debug)]
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
            intercept: self.intercept,
        }
    }
}

#[derive(Clone)]
#[pyclass]
pub struct PyLinearModelParams {
    #[pyo3(get, set)]
    pub gradient: f64,
    #[pyo3(get, set)]
    pub intercept: f64,
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
            corr: self.corr,
        }
    }
}

#[derive(Clone)]
#[pyclass]
pub struct PyLinearModelEval {
    #[pyo3(get, set)]
    pub error: f64,
    #[pyo3(get, set)]
    pub corr: f64,
}


pub struct LinearModel {
    pub params: LinearModelParams,
    pub eval: LinearModelEval,
}

impl LinearModel {
    pub fn fit_from_params(
        params: LinearModelParams,
        model_error: LinearModelError,
        model_corr: LinearModelCorr,
        x: &[f64],
        y: &[f64],
    ) -> Self {

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
            eval: self.eval.to_python(),
        }
    }
}

#[derive(Clone)]
#[pyclass]
pub struct PyLinearModel {
    #[pyo3(get, set)]
    pub params: PyLinearModelParams,
    #[pyo3(get, set)]
    pub eval: PyLinearModelEval,
}


pub struct LinearModelNew<'a> {
    x: &'a [f64],
    y: &'a [f64],

    // The slopes for each j that can be computed
    pub slopes_j: Vec<Vec<f64>>,
    pub slopes_nonzero_mask: Vec<BitVec>,
    pub slopes_j_sorted: Vec<Vec<usize>>,
}


impl<'a> LinearModelNew<'a> {
    /// Create a new Linear Model and setup some values that we
    /// don't want to compute again after jackknifing
    pub fn new(x: &'a [f64], y: &'a [f64]) -> DeltaPDResult<Self> {

        // Ensure the vectors are the same length
        if x.len() != y.len() {
            return Err(DeltaPDError::Error("X and Y vectors are not the same length.".to_string()));
        }

        // Compute values that will not change
        let n = x.len();
        let mut slopes_nonzero_mask: Vec<BitVec> = Vec::with_capacity(n);
        let mut slopes_j: Vec<Vec<f64>> = Vec::with_capacity(n);
        let mut slopes_j_sorted: Vec<Vec<usize>> = Vec::with_capacity(n);

        // Perform the first part of the repeated-median calculation
        for i in 0..n {
            let mut cur_slopes_nonzero_mask = bitvec![0; n];
            let mut cur_slopes_j: Vec<f64> = vec![0.0; n];
            for j in 0..n {
                let dx = x[i] - x[j];
                let dy = y[i] - y[j];

                // Skip data points that will have a zero denominator
                if dx != 0.0 {
                    cur_slopes_nonzero_mask.set(j, true);
                    cur_slopes_j[j] = dy / dx;
                }
            }

            // Sort the slopes and get the indices
            let cur_slopes_j_sorted = argsort_by_vec(&cur_slopes_j);

            // Save the values for this model
            slopes_nonzero_mask.push(cur_slopes_nonzero_mask);
            slopes_j.push(cur_slopes_j);
            slopes_j_sorted.push(cur_slopes_j_sorted);
        }

        Ok(Self {
            x,
            y,
            slopes_j,
            slopes_nonzero_mask,
            slopes_j_sorted,
        })
    }

    /// Calculate the remainder of the repeated median method.
    /// Optionally exclude datapoints indicated by the presence_bv.
    pub fn compute(&self, presence_bv: Option<&BitVec>) -> LinearModelParams {
        let n = self.x.len();

        let mut median_values: Vec<f64> = Vec::with_capacity(n);

        // Go over each row
        for i in 0..n {
            let mut sorted_values: Vec<f64> = Vec::with_capacity(n);

            // Get the current mask for the row
            let cur_row_mask: BitVec = {
                // Get the current row bitmask
                let cur_row_bitmask = &self.slopes_nonzero_mask[i];

                // If a user has specified a taxon then do a bitwise comparison
                if let Some(cur_label_bitmask) = presence_bv {
                    bitvec_boolean_a_and_not_b(cur_row_bitmask, cur_label_bitmask)
                } else {
                    // TODO: Remove this clone if possible?
                    cur_row_bitmask.clone()
                }
            };

            // Go over each column in the sorted order
            for j in 0..n {
                let cur_idx = self.slopes_j_sorted[i][j];
                let cur_mask = *cur_row_mask.get(cur_idx).as_deref().unwrap();
                if cur_mask {
                    let cur_value = self.slopes_j[i][cur_idx];
                    sorted_values.push(cur_value);
                }
            }

            // Calculate the median
            let median = calc_median_sorted(&sorted_values);
            median_values.push(median);
        }

        let medslope = calc_median(&median_values);

        let mut medinter_vec: Vec<f64> = Vec::with_capacity(n);

        if let Some(cur_bv) = presence_bv {
            // If filtering is applied, then get only those points that do not contain the taxon
            for i in cur_bv.iter_zeros() {
                medinter_vec.push(self.y[i] - medslope * self.x[i]);
            }
        } else {
            // Otherwise, compute across the whole dataset
            for i in 0..n {
                medinter_vec.push(self.y[i] - medslope * self.x[i]);
            }
        }

        let medinter = calc_median(&medinter_vec);
        LinearModelParams::new(medslope, medinter)
    }

    /// Fit the model and return the error and correlation
    /// Optionally exclude datapoints indicated by the presence_bv.
    pub fn fit(&self, presence_bv: Option<&BitVec>, params: LinearModelParams, model_error: LinearModelError, model_corr: LinearModelCorr) -> LinearModelEval {
        let error_fn = model_error.get_fn();
        let corr_fn = model_corr.get_fn();

        // If the BitVec has been provided, subset the x/y data to those points
        if let Some(taxon_bv) = presence_bv {
            let mut y_hat: Vec<f64> = Vec::with_capacity(taxon_bv.len());
            let mut y: Vec<f64> = Vec::with_capacity(taxon_bv.len());
            for i in taxon_bv.iter_zeros() {
                y_hat.push(params.gradient * self.x[i] + params.intercept);
                y.push(self.y[i]);
            }
            let error = error_fn(&y, &y_hat);
            let corr = corr_fn(&y, &y_hat);
            LinearModelEval::new(error, corr)
        } else {
            let mut y_hat = Vec::with_capacity(self.x.len());
            for i in 0..self.x.len() {
                y_hat.push(params.gradient * self.x[i] + params.intercept);
            }
            let error = error_fn(self.y, &y_hat);
            let corr = corr_fn(self.y, &y_hat);
            LinearModelEval::new(error, corr)
        }
    }
}


#[test]
fn test_linear_model_new() {
    use crate::util::numeric::round_f64;

    let x = vec![0.09, -0.4, -0.03, -0.33, -0.1, -0.48, 0.12, -0.01, 0.49, -0.47];
    let y = vec![0.09, -0.2, -0.05, -0.3, -0.45, -0.08, -0.39, -0.42, 0.01, -0.19];

    let model = LinearModelNew::new(&x, &y).unwrap();
    let results = model.compute(None);

    assert_eq!(round_f64(results.gradient, 5), -0.33898);
    assert_eq!(round_f64(results.intercept, 5), -0.34246);

    let bit_vec = bitvec![0, 0, 1, 0, 0, 1, 0, 0, 0, 0];
    let results = model.compute(Some(&bit_vec));
    assert_eq!(round_f64(results.gradient, 5), 0.01538);
    assert_eq!(round_f64(results.intercept, 5), -0.24438);
}
