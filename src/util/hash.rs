use std::fs::File;
use std::io::{BufReader, Read};
use std::path::PathBuf;

use md5::{Digest, Md5};
use pyo3::exceptions::PyValueError;
use pyo3::{pyfunction, PyResult};

use crate::model::error::{DeltaPDError, DeltaPDResult};

pub fn file_md5(path: &PathBuf) -> DeltaPDResult<String> {
    let file = File::open(path).map_err(|e| DeltaPDError::IoError(e))?;
    let mut reader = BufReader::new(file);
    let mut hasher = Md5::new();
    let mut buffer = [0; 1024];

    loop {
        let count = reader.read(&mut buffer).map_err(|e| DeltaPDError::IoError(e))?;
        if count == 0 { break; }
        hasher.update(&buffer[..count]);
    }

    let result = hasher.finalize();
    let hex = format!("{:x}", result);
    Ok(hex)
}

#[pyfunction]
#[pyo3(name = "file_md5")]
pub fn file_md5_py(path: String) -> PyResult<String> {
    let hash = file_md5(&PathBuf::from(&path));
    match hash {
        Ok(h) => Ok(h),
        Err(e) => Err(PyValueError::new_err(format!("Unable to calculate MD5 {e}: {}", &path)))
    }
}

#[test]
fn test_file_md5() {
    use tempfile::tempdir;
    use std::io::Write;

    // Create a temporary directory
    let dir = tempdir().expect("Failed to create temp dir");
    let file_path = dir.path().join("test.txt");
    let mut file = File::create(&file_path).expect("Unable to create file");
    writeln!(file, "ae9fja0921").expect("Unable to write to file");

    // Calculate expected MD5 hash of "foo"
    let expected_md5 = "dfb19f01c983e75abbc842eeef38697e";

    // Call the function and compare
    let result = file_md5(&file_path).expect("Failed to get MD5 hash");
    assert_eq!(result, expected_md5);

    // Clean up
    dir.close().expect("Failed to clean up temp dir");
}