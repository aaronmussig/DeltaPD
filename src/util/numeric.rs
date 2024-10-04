pub fn round_f64(x: f64, decimals: u32) -> f64 {
    let factor = 10_f64.powi(decimals as i32);
    (x * factor).round() / factor
}

#[test]
fn test_round() {
    assert_eq!(round_f64(1.23456789, 0), 1.0);
    assert_eq!(round_f64(1.23456789, 1), 1.2);
    assert_eq!(round_f64(1.23456789, 3), 1.235);
    assert_eq!(round_f64(1.23456789, 4), 1.2346);
    assert_eq!(round_f64(1.23456789, 5), 1.23457);
    assert_eq!(round_f64(1.23456789, 6), 1.234568);
    assert_eq!(round_f64(1.9, 0), 2.0);
}

