use ndarray::{Array1, ArrayView1};

// pub fn get_nonzero(x: &ArrayView1<f64>) -> Array1<f64> {
//     let mut nonzero: Vec<f64> = Vec::new();
//     for i in 0..x.len() {
//         if x[i] != 0.0 {
//             nonzero.push(x[i]);
//         }
//     }
//     return Array1::from_vec(nonzero);
// }

pub fn apply_mask_1d(x: &ArrayView1<f64>, mask: &ArrayView1<bool>) -> Array1<f64> {
    let mut masked: Vec<f64> = Vec::with_capacity(x.len());
    for i in 0..x.len() {
        if mask[i] {
            masked.push(x[i]);
        }
    }
    return Array1::from_vec(masked);
}