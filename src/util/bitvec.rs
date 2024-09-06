use bitvec::prelude::*;


/// For two BitVecs, compute (a & !b)
pub fn bitvec_boolean_a_and_not_b(a: &BitVec, b: &BitVec) -> BitVec {
    assert_eq!(a.len(), b.len(), "BitVecs must have the same length");

    let a_slice = a.as_raw_slice();
    let b_slice = b.as_raw_slice();

    // Perform a AND (NOT b) on each usize element
    let mut result_storage = vec![0; a_slice.len()];
    for (i, (a_val, b_val)) in a_slice.iter().zip(b_slice.iter()).enumerate() {
        result_storage[i] = a_val & !b_val; // AND a with the negation of b
    }

    let mut result = BitVec::from_vec(result_storage);
    result.truncate(a.len());
    result
}

#[test]
fn test_bitvec_boolean_a_and_not_b() {
    use rand::{thread_rng, Rng};

    let n_bits = 10_000;
    let mut a_vec = bitvec![0; n_bits];
    let mut b_vec = bitvec![0; n_bits];
    let mut expected_bv = bitvec![0; n_bits];

    let mut rng = thread_rng();

    for i in 0..n_bits {
        let a_val: bool = rng.gen();
        let b_val: bool = rng.gen();
        let expected_val = a_val & !b_val;
        a_vec.set(i, a_val);
        b_vec.set(i, b_val);
        expected_bv.set(i, expected_val);
    }

    let result = bitvec_boolean_a_and_not_b(&a_vec, &b_vec);
    assert_eq!(result, expected_bv);
}