use ndarray::Array2;
use crate::model::error::{DeltaPDError, DeltaPDResult};

pub fn contains_nan(x: &Array2<f64>) -> bool {
    for v in x.iter() {
        if v.is_nan() {
            return true;
        }
    }
    false
}

pub fn exit_if_nan(x: &Array2<f64>) -> DeltaPDResult<()> {
    if contains_nan(x) {
        return Err(DeltaPDError::Error("NaN found in matrix".to_string()));
    }
    Ok(())
}