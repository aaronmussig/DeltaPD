use std::cmp::Ordering;

use ndarray::Data;
use ndarray::{s, ArrayBase, ArrayView, Ix1, Ix2};
use rayon::prelude::*;

/// https://github.com/rust-ndarray/ndarray/issues/1145
pub fn argsort_by<S, F>(arr: &ArrayBase<S, Ix1>, mut compare: F) -> Vec<usize>
where
    S: Data,
    F: FnMut(&S::Elem, &S::Elem) -> Ordering,
{
    let mut indices: Vec<usize> = (0..arr.len()).collect();
    indices.sort_unstable_by(move |&i, &j| compare(&arr[i], &arr[j]));
    indices
}

/// Returns the indices that would sort the array.
pub fn argsort_by_vec(arr: &[f64]) -> Vec<usize> {
    let mut indices: Vec<usize> = (0..arr.len()).collect();
    indices.sort_unstable_by(|&i, &j| arr[i].partial_cmp(&arr[j]).expect("Elements must not be NaN."));
    indices
}

#[test]
fn test_argsort_by_vec() {
    let arr = vec![1.0, 3.0, 2.0];
    let indices = argsort_by_vec(&arr);
    assert_eq!(indices, vec![0, 2, 1]);

    let arr = vec![1.0, 3.0, 2.0, 4.0];
    let indices = argsort_by_vec(&arr);
    assert_eq!(indices, vec![0, 2, 1, 3]);
}

/// Returns the nearest neighbours for each taxon in a distance matrix.
pub fn get_nn_from_distance_matrix(matrix: &ArrayView<f64, Ix2>) -> Vec<Vec<usize>> {
    let n_slices = matrix.shape()[0];
    let knn_idx = (0..n_slices).into_par_iter().map(|n| {
        let slice = matrix.slice(s![n, ..]);
        argsort_by(&slice, |a, b| {
            a.partial_cmp(b).expect("Elements must not be NaN.")
        })
    }).collect::<Vec<Vec<usize>>>();
    knn_idx
}


// TODO: Write unit tests for get_nn_from_distance_matrix
