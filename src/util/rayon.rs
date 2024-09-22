use crate::model::error::DeltaPDError;

pub fn create_pool(num_threads: usize) -> Result<rayon::ThreadPool, DeltaPDError> {
    match rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
    {
        Err(e) => Err(DeltaPDError::RayonError(e)),
        Ok(pool) => Ok(pool),
    }
}