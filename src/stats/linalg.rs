use ndarray::Array2;

use crate::model::linalg::LinearModelParams;
use crate::ndarray::filter::apply_mask_1d;
use crate::stats::vec::{calc_mean, calc_median, calc_stddev};

/// Calculate the gradient of X/Y coordinates using the Theil-Sen method.
/// Assumes the intercept is 0.
pub fn calc_theil_sen_gradient(x: &[f64], y: &[f64]) -> LinearModelParams {
    assert_eq!(
        x.len(),
        y.len(),
        "Theil-Sen gradient requires equal length vectors."
    );

    // Compute the estimate at each data point
    let mut estimate = Vec::with_capacity(x.len());
    for i in 0..x.len() {
        let x_i = x[i];
        let y_i = y[i];
        if x_i != 0.0 {
            estimate.push(y_i / x_i);
        }
    }

    // Take the median value
    LinearModelParams::new(calc_median(&estimate), 0.0)
}


// TODO: Implement O(n) algorithm: 10.1016/S0020-0190(03)00350-8
pub fn calc_repeated_median(x: &[f64], y: &[f64]) -> LinearModelParams {
    let n = x.len();
    let mut deltax: Array2<f64> = Array2::zeros((n, n));
    let mut deltay: Array2<f64> = Array2::zeros((n, n));

    for i in 0..n {
        for j in 0..n {
            deltax[[i, j]] = x[i] - x[j];
            deltay[[i, j]] = y[i] - y[j];
        }
    }

    let mut slopes: Vec<f64> = Vec::with_capacity(n);

    for j in 0..n {
        let id_nonzero = deltax.row(j).mapv(|x| x != 0.0);
        let deltax_nonzero = apply_mask_1d(&deltax.row(j), &id_nonzero.view());
        let deltay_nonzero = apply_mask_1d(&deltay.row(j), &id_nonzero.view());

        let slopes_j: Vec<f64> = (&deltay_nonzero / &deltax_nonzero).to_vec();
        let medslope_j = calc_median(&slopes_j);
        slopes.push(medslope_j);
    }

    let medslope = calc_median(&slopes);

    let medinter_vec: Vec<f64> = x
        .iter()
        .zip(y.iter())
        .map(|(&x_i, &y_i)| y_i - medslope * x_i)
        .collect();
    let medinter = calc_median(&medinter_vec);

    LinearModelParams::new(medslope, medinter)
}

#[test]
fn test_repeated_median() {
    use crate::util::numeric::round_f64;

    let x = vec![0.09, -0.4, -0.03, -0.33, -0.1, -0.48, 0.12, -0.01, 0.49, -0.47];
    let y = vec![0.09, -0.2, -0.05, -0.3, -0.45, -0.08, -0.39, -0.42, 0.01, -0.19];

    let results = calc_repeated_median(&x, &y);
    assert_eq!(round_f64(results.gradient, 5), -0.33898);
    assert_eq!(round_f64(results.intercept, 5), -0.34246);
}


/// Calculate the coefficient of determination
pub fn calc_r2(y: &[f64], y_hat: &[f64]) -> f64 {
    assert_eq!(y.len(), y_hat.len(), "R2 requires equal length vectors.");

    let y_mean = calc_mean(y);

    let (numerator, denominator) = y.iter().zip(y_hat.iter()).fold(
        (0.0, 0.0),
        |(num, denom), (&y_i, &y_hat_i)| {
            (
                num + (y_i - y_hat_i).powi(2),
                denom + (y_i - y_mean).powi(2),
            )
        },
    );

    1.0 - (numerator / denominator)
}

#[test]
fn test_calc_r2() {
    assert_eq!(
        calc_r2(&[84.0, 0.0, 62.0, 20.0], &[0.0, 59.04, 11.48, 137.76]),
        -5.112312310133756
    );
    assert_eq!(calc_r2(&[22.0, 0.0], &[0.0, 103.7]), -45.43673553719008);
}

/// Calculate the mean squared error.
pub fn calc_mse(y: &[f64], y_hat: &[f64]) -> f64 {
    assert_eq!(y.len(), y_hat.len(), "MSE requires equal length vectors.");

    y.iter()
        .zip(y_hat.iter())
        .map(|(&y_i, &y_hat_i)| (y_i - y_hat_i).powi(2))
        .sum::<f64>()
        / y.len() as f64
}

#[test]
fn test_calc_mse() {
    assert_eq!(
        calc_mse(&[58.0, 67.0, 18.0], &[35.96, 0.0, 22.62]),
        1665.3686666666665
    );
}

/// Calculate the root mean squared error.
pub fn calc_rmse(y: &[f64], y_hat: &[f64]) -> f64 {
    calc_mse(y, y_hat).sqrt()
}

#[test]
fn test_calc_rmse() {
    use crate::util::numeric::round_f64;

    assert_eq!(
        round_f64(calc_rmse(&[58.0, 67.0, 18.0], &[35.96, 0.0, 22.62]), 5),
        40.80893
    );
}

/// Calculate the mean squared error normalised by the variance.
pub fn calc_mse_norm(y: &[f64], y_hat: &[f64]) -> f64 {
    // Here we can call RMSE instead of MSE then square rooting as it's the same
    calc_rmse(y, y_hat) / calc_stddev(y)
}

#[test]
fn test_calc_mse_norm() {
    assert_eq!(
        calc_mse_norm(&[1.0, 38.0, 44.0, 57.0], &[52.8, 91.2, 137.6, 9.6]),
        3.0902793708005385
    );
}

// Calculate the pearson correlation coefficient
pub fn calc_pearson_correlation(x: &[f64], y: &[f64]) -> f64 {
    assert_eq!(
        x.len(),
        y.len(),
        "Pearson correlation requires equal length vectors."
    );

    let x_mean = calc_mean(x);
    let y_mean = calc_mean(y);
    let (mut numerator, mut denominator_x, mut denominator_y) = (0.0, 0.0, 0.0);

    for i in 0..x.len() {
        let (x_diff, y_diff) = (x[i] - x_mean, y[i] - y_mean);
        numerator += x_diff * y_diff;
        denominator_x += x_diff.powi(2);
        denominator_y += y_diff.powi(2);
    }

    numerator / (denominator_x * denominator_y).sqrt()
}


#[test]
fn test_calc_calc_pearson_correlation() {
    use crate::util::numeric::round_f64;

    assert_eq!(
        round_f64(
            calc_pearson_correlation(&[5.0, 4.0, 5.0, 2.0, 4.0, ], &[0.0, 6.0, 3.0, 7.0, 0.0, ]),
            6,
        ),
        round_f64(-0.6864282924115621, 6)
    );
    assert_eq!(
        calc_pearson_correlation(&[-2.0, -1.0, 0.0, 1.0, 2.0], &[-2.0, -1.0, 0.0, 1.0, 2.0]),
        1.0
    );
}

/// Compute the estimate of a vector of X values.
pub fn calc_y_hat(x: &[f64], m: f64) -> Vec<f64> {
    x.iter().map(|&x_i| m * x_i).collect()
}


#[test]
fn test_calc_y_hat() {
    assert_eq!(calc_y_hat(&[1.0, 2.0, 3.0], 2.0), vec![2.0, 4.0, 6.0]);
}
