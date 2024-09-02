use std::path::PathBuf;

use crate::model::linalg::LinearModelParams;

pub type ModelFn = fn(&[f64], &[f64]) -> LinearModelParams;


pub type CorrFn = fn(&[f64], &[f64]) -> f64;

pub type ErrorFn = fn(&[f64], &[f64]) -> f64;


#[derive(Debug, Eq, Hash, PartialEq)]
pub struct SequenceId {
    pub value: String,
}

#[derive(Debug)]
pub struct GenomeId {
    pub value: String,
}

#[derive(Debug)]
pub struct ContigLength {
    pub value: u32,
}


#[derive(Debug)]
pub struct Knn {
    pub value: usize,
}

#[derive(Debug)]
pub struct MaxThreads {
    pub value: usize,
}

#[derive(Debug)]
pub struct PathRefTree {
    pub value: PathBuf,
}

#[derive(Debug)]
pub struct PathQryTree {
    pub value: PathBuf,
}
