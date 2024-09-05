use phylodm::tree::Taxon;

pub struct Job<'a>  {
    pub taxon: &'a Taxon,
    pub replicate: usize,
}

impl<'a> Job<'a> {
    pub fn new(taxon: &'a Taxon, replicate: usize) -> Self {
        Job {
            taxon,
            replicate,
        }
    }
}