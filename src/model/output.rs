use phylodm::tree::Taxon;

pub struct DeltaPdOutput {
    pub taxon: Taxon,
    pub replicate: usize,
    pub rel_influence: f64,
    pub prop_std_error: f64,
    pub gradient: f64,
    pub intercept: f64,
}


impl DeltaPdOutput {
    pub fn new(taxon: &Taxon, replicate: usize, rel_influence: f64, prop_std_error: f64, gradient: f64, intercept: f64) -> Self {
        Self {
            taxon: taxon.clone(),
            replicate,
            rel_influence,
            prop_std_error,
            gradient,
            intercept,
        }
    }
}


