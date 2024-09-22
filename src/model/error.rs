
#[derive(Debug)]
pub enum DeltaPDError {
    Error(String),
    Exit(),
    CsvError(csv::Error),
    IoError(std::io::Error),
    RayonError(rayon::ThreadPoolBuildError),
    PhyloDMError(phylodm::error::PhyloErr),
    ParseError(std::string::FromUtf8Error),
    ParseFloatError(std::num::ParseFloatError),
    ParseIntError(std::num::ParseIntError),
    ThreadPoolBuildError(rayon::ThreadPoolBuildError),
}

impl std::fmt::Display for DeltaPDError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            DeltaPDError::Error(e) => write!(f, "{}", e),
            DeltaPDError::CsvError(e) => write!(f, "CSV error: {}", e),
            DeltaPDError::Exit() => write!(f, "Exit"),
            DeltaPDError::IoError(e) => write!(f, "IO error: {}", e),
            DeltaPDError::RayonError(e) => write!(f, "Rayon error: {}", e),
            DeltaPDError::PhyloDMError(e) => write!(f, "PhyloDM error: {}", e),
            DeltaPDError::ParseError(e) => write!(f, "Parse error: {}", e),
            DeltaPDError::ParseFloatError(e) => write!(f, "Parse float error: {}", e),
            DeltaPDError::ParseIntError(e) => write!(f, "Parse int error: {}", e),
            DeltaPDError::ThreadPoolBuildError(e) => write!(f, "ThreadPoolBuildError: {}", e),
        }
    }
}

impl std::error::Error for DeltaPDError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            DeltaPDError::Error(_) => None,
            DeltaPDError::Exit() => None,
            DeltaPDError::CsvError(e) => Some(e),
            DeltaPDError::IoError(e) => Some(e),
            DeltaPDError::RayonError(e) => Some(e),
            DeltaPDError::PhyloDMError(e) => Some(e),
            DeltaPDError::ParseError(e) => Some(e),
            DeltaPDError::ParseFloatError(e) => Some(e),
            DeltaPDError::ParseIntError(e) => Some(e),
            DeltaPDError::ThreadPoolBuildError(e) => Some(e),
        }
    }
}


pub type DeltaPDResult<T> = Result<T, DeltaPDError>;