use std::path::PathBuf;

use crate::method::deltapd::run_deltapd;
use crate::model::metadata::MetadataFile;
use crate::model::params::{Params, PyParams};
use crate::model::pdm::{PyDistMatrix, QryDistMatrix, RefDistMatrix};
use crate::model::result::{OutputResultSmall, PyOutputResultSmall};
use pyo3::prelude::*;

pub struct DeltaPD {
    pub qry_dm: QryDistMatrix,
    pub ref_dm: RefDistMatrix,
    pub metadata: MetadataFile,
}

impl DeltaPD {
    pub fn new(qry_dm: QryDistMatrix, ref_dm: RefDistMatrix, metadata: MetadataFile) -> Self {
        Self {
            qry_dm,
            ref_dm,
            metadata,
        }
    }

    pub fn run(&self, params: &Params) -> Vec<OutputResultSmall> {
        run_deltapd(
            &self.qry_dm,
            &self.ref_dm,
            &self.metadata,
            &params,
        ).unwrap()
    }
}

#[pyclass]
pub struct PyDeltaPD {
    pub deltapd: DeltaPD,
}

#[pymethods]
impl PyDeltaPD {
    #[new]
    pub fn new(qry_dm: PyDistMatrix, ref_dm: PyDistMatrix, metadata_path: &str, delimiter: &str) -> Self {
        let metadata_path = PathBuf::from(metadata_path);
        let delim = delimiter.as_bytes()[0];
        let metadata_file = MetadataFile::read(&metadata_path, delim).expect("TODO");
        Self {
            deltapd: DeltaPD::new(
                QryDistMatrix { dm: qry_dm.dm },
                RefDistMatrix { dm: ref_dm.dm },
                metadata_file,
            )
        }
    }

    pub fn run(&self, params: &PyParams) -> PyResult<Vec<PyOutputResultSmall>> {
        let result = self.deltapd.run(&params.params);
        let mut output = Vec::with_capacity(result.len());
        for res in result {
            output.push(res.to_python());
        }
        Ok(output)
    }
}


