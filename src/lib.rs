use pyo3::{Bound, pymodule, PyResult, wrap_pyfunction};
use pyo3::prelude::{PyModule, PyModuleMethods};
use crate::model::linalg::{PyLinearModelCorr, PyLinearModelError, PyLinearModelType};
use crate::model::params::PyParams;
use crate::model::pdm::PyDistMatrix;
use crate::python::deltapd::PyDeltaPD;
use crate::util::hash::file_md5_py;

pub mod method;
pub mod ndarray;

pub mod python;

pub mod util;
pub mod model;
pub mod stats;



/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn rs_deltapd(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyDeltaPD>()?;
    m.add_class::<PyDistMatrix>()?;
    m.add_class::<PyParams>()?;
    m.add_class::<PyLinearModelType>()?;
    m.add_class::<PyLinearModelError>()?;
    m.add_class::<PyLinearModelCorr>()?;
    m.add_function(wrap_pyfunction!(file_md5_py, m)?)?;
    Ok(())
}
