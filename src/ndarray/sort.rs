use std::cmp::Ordering;

use ndarray::{ArrayBase, ArrayView, Ix1, Ix2, s};
use ndarray::Data;
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
