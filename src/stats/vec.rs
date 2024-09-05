/// Calculate the median value of a vector.
/// // TODO: Can make this mutable for a bit more performance gain
pub fn calc_median(x: &[f64]) -> f64 {
    if x.is_empty() {
        panic!("Cannot calculate the median of an empty vector.");
    }
    let mut sorted = x.to_vec();
    sorted.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let mid = sorted.len() / 2;

    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) / 2.0
    } else {
        sorted[mid]
    }
}

pub fn calc_median_sorted(x: &[f64]) -> f64 {
    if x.is_empty() {
        panic!("Cannot calculate the median of an empty vector.");
    }

    let n = x.len();
    let mid = n / 2;

    if n % 2 == 0 {
        (x[mid - 1] + x[mid]) / 2.0
    } else {
        x[mid]
    }
}


#[test]
fn test_calc_median() {
    assert_eq!(calc_median(&[5.0, 3.0, 4.0, 2.0, 1.0]), 3.0);
    assert_eq!(calc_median(&[1.0, 3.0, 4.0, 2.0]), 2.5);
}

/// Calculate the mean of a slice.
pub fn calc_mean(x: &[f64]) -> f64 {
    let total: f64 = x.iter().sum();
    total / (x.len() as f64)
}

#[test]
fn test_calc_mean() {
    assert_eq!(calc_mean(&[1.0, 2.0, 3.0]), 2.0);
    assert_eq!(calc_mean(&[5.0, 5.0]), 5.0);
    assert_eq!(calc_mean(&[10.0, 20.0]), 15.0);
}

/// Calculate the standard deviation of a slice.
pub fn calc_stddev(x: &[f64]) -> f64 {
    let mean = calc_mean(x);
    let total: f64 = x.iter().map(|&x_i| (x_i - mean).powi(2)).sum();
    (total / x.len() as f64).sqrt()
}

#[test]
fn test_calc_stddev() {
    assert_eq!(calc_stddev(&[1.0, 2.0, 3.0, 4.0]), 1.118033988749895);
    assert_eq!(calc_stddev(&[1.0, 2.0]), 0.5);
}
