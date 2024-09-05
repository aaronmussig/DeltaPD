use std::collections::{HashMap, HashSet};
use crate::model::error::{DeltaPDError, DeltaPDResult};
use crate::model::linalg::{LinearModel, LinearModelCorr, LinearModelError, LinearModelType};
use crate::stats::linalg::{calc_rmse, calc_theil_sen_gradient, calc_y_hat};
use crate::stats::vec::calc_mean;

/// Determine the relative influence of each data point in the jackknife set. Efron (1992).
/// Equation 2.4 of https://www.jstor.org/stable/2345949
pub fn calc_relative_jackknife_influence(x: &[f64]) -> Vec<f64> {
    let x_bar = calc_mean(x);
    let n = x.len() as f64;

    let sum_squares: f64 = x.iter()
        .map(|&x_i| ((n - 1.0) * (x_bar - x_i)).powi(2))
        .sum();
    let denominator = (sum_squares / (n - 1.0)).sqrt();

    x.iter()
        .map(|&x_i| (n - 1.0) * (x_bar - x_i) / denominator)
        .collect()
}


/// Calculate the proportion of the std error each point contributes.
pub fn calc_prop_std_error_of_jackknife(x: &[f64]) -> Vec<f64> {
    let n_inv = 1.0 / (x.len() as f64 - 1.0);
    x.iter()
        .map(|&x_i| x_i.powi(2) * n_inv)
        .collect()
}

/// Jackknife a set of data using the specified function.
pub fn jackknife(x: &[f64], y: &[f64], f: fn(&[f64], &[f64]) -> f64) -> Vec<f64> {
    assert_eq!(x.len(), y.len(), "Jacknife requires equal length vectors.");
    let n = x.len();
    let mut results = Vec::with_capacity(n);
    for i in 0..n {
        let (x_left, x_right) = x.split_at(i);
        let (y_left, y_right) = y.split_at(i);
        let x_copy = [x_left, &x_right[1..]].concat();
        let y_copy = [y_left, &y_right[1..]].concat();
        results.push(f(&x_copy, &y_copy));
    }
    results
}

#[test]
fn test_jackknife() {
    fn _sum_two_vecs(x: &[f64], y: &[f64]) -> f64 {
        let mut tot = 0.0;
        for i in 0..x.len() {
            tot += x[i] + y[i];
        }
        tot
    }
    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let results = jackknife(&x, &y, _sum_two_vecs);
    assert_eq!(results, vec![28.0, 26.0, 24.0, 22.0, 20.0]);
}

#[test]
fn test_relative_jackknife_influence() {
    use crate::stats::linalg::calc_pearson_correlation;
    use crate::util::numeric::round_f64;

    // law school data set points efron 1983 - american statistician vol 37 no 1.
    let law_z = [
        3.39, 3.30, 2.81, 3.03, 3.44, 3.07, 3.00, 3.43, 3.36, 3.13, 3.12, 2.74, 2.76, 2.88, 2.96,
    ];
    let law_y = [
        576.0, 635.0, 558.0, 578.0, 666.0, 580.0, 555.0, 661.0, 651.0, 605.0, 653.0, 575.0, 545.0,
        572.0, 594.0,
    ];

    // Pearson correlation
    let law_jackknifed_data = jackknife(&law_z, &law_y, calc_pearson_correlation);
    let law_jackknife_expected = [
        -2.97, 0.31, 0.53, 0.0, 1.13, -0.1, -0.22, 1.01, 0.61, -0.01, -1.07, -0.25, 0.9, 0.22, -0.1,
    ];
    let mut law_jackknife_influence = vec![0.0; law_jackknifed_data.len()];
    for (i, u_i) in calc_relative_jackknife_influence(&law_jackknifed_data)
        .into_iter()
        .enumerate()
    {
        law_jackknife_influence[i] = round_f64(u_i, 2);
    }
    assert_eq!(law_jackknife_influence, law_jackknife_expected);

    let prop_contrib = calc_prop_std_error_of_jackknife(&law_jackknife_expected);
    assert_eq!(round_f64(prop_contrib[0], 2), 0.63);

    let mut prop_sum = 0.0;
    for prop_i in prop_contrib {
        prop_sum += prop_i
    }
    assert_eq!(round_f64(prop_sum, 2), 1.0);
}

pub struct JackknifeResult {
    pub idx: usize,
    pub grad: f64,
    pub rmse: f64,
}

pub fn jackknife_theil_sen_rmse(x: &[f64], y: &[f64]) -> Vec<JackknifeResult> {
    assert_eq!(x.len(), y.len(), "Jacknife requires equal length vectors.");
    let mut results: Vec<JackknifeResult> = Vec::with_capacity(x.len());

    for i in 0..x.len() {
        let mut x_copy = x.to_vec();
        x_copy.remove(i);

        let mut y_copy = y.to_vec();
        y_copy.remove(i);

        let lin_model = calc_theil_sen_gradient(&x_copy, &y_copy);
        let grad = lin_model.gradient;
        let y_hat = calc_y_hat(&x_copy, grad);
        let rmse = calc_rmse(&y_copy, &y_hat);

        results.push(JackknifeResult { idx: i, grad, rmse });
    }
    results
}

pub fn jackknife_fn(
    x: &[f64],
    y: &[f64],
    model_type: LinearModelType,
    model_error: LinearModelError,
    model_corr: LinearModelCorr,
) -> DeltaPDResult<Vec<LinearModel>> {
    // Validate the input
    if x.len() != y.len() {
        return Err(DeltaPDError::Error(
            "Jackknife requires equal length vectors.".to_string(),
        ));
    }

    let n = x.len();
    let mut results = Vec::with_capacity(n);

    // Where x is the reference data points and y are the query data points

    for i in 0..n {
        let (x_left, x_right) = x.split_at(i);
        let (y_left, y_right) = y.split_at(i);
        let x_copy = [x_left, &x_right[1..]].concat();
        let y_copy = [y_left, &y_right[1..]].concat();

        let linear_model = LinearModel::fit(model_type, model_error, model_corr, &x_copy, &y_copy);

        results.push(linear_model);
    }
    Ok(results)
}


pub fn jackknife_fn_masked(
    x: &[f64],
    y: &[f64],
    model_type: LinearModelType,
    model_error: LinearModelError,
    model_corr: LinearModelCorr,
    mask: &Vec<(usize, usize)>
) -> DeltaPDResult<(Vec<LinearModel>, Vec<usize>)> {

    // Validate the input
    if x.len() != y.len() {
        return Err(DeltaPDError::Error(
            "Jackknife requires equal length vectors.".to_string(),
        ));
    }

    // Create a mapping of the taxon ID to the indices
    let mut taxon_to_indices: HashMap<usize, HashSet<usize>> = HashMap::new();
    for (i, (x_i, y_i)) in mask.iter().enumerate() {
        taxon_to_indices.entry(*x_i).or_insert(HashSet::new()).insert(i);
        taxon_to_indices.entry(*y_i).or_insert(HashSet::new()).insert(i);
    }

    let n = taxon_to_indices.len();
    let mut results = Vec::with_capacity(n);

    let mut taxon_out = Vec::new();
    for (taxon, indices) in taxon_to_indices.iter() {
        taxon_out.push(*taxon);
        let mut new_x: Vec<f64> = Vec::with_capacity(x.len()-indices.len());
        let mut new_y: Vec<f64> = Vec::with_capacity(x.len()-indices.len());

        // get all points without this taxon
        for i in 0..x.len() {
            if !indices.contains(&i) {
                new_x.push(x[i]);
                new_y.push(y[i]);
            }
        }

        let linear_model = LinearModel::fit(model_type, model_error, model_corr, &new_x, &new_y);
        results.push(linear_model);
    }

    Ok((results, taxon_out))
}