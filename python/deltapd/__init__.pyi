from enum import Enum
from typing import Set

from deltapd.model.params import Direction


class PyOutputResultSmall:
    error_mean: float
    error_median: float
    error_std: float
    taxon: str


class PyDeltaPD:

    def __init__(self, qry_dm: PyDistMatrix, ref_dm: PyDistMatrix, metadata_path: str, delimiter: str):
        pass

    def run(self, params: PyParams) -> list[PyOutputResultSmall]:
        pass


class PyDistMatrix:

    def __init__(self, taxa: tuple[str, int], edges: tuple[int, int, float]):
        pass


def file_md5(path: str) -> str:
    pass


class PyLinearModelType(Enum):
    RepeatedMedian = 0
    TheilSen = 1


class PyLinearModelError(Enum):
    MSE = 0
    NormMSE = 1
    RMSE = 2


class PyLinearModelCorr(Enum):
    R2 = 0
    Pearson = 1


class PyParams:
    def __init__(self, cpus: int, direction: Direction, sample_size: float, knn: int, replicates: int, taxa: Set[str], model: PyLinearModelType,
                 error: PyLinearModelError, corr: PyLinearModelCorr, debug: bool, output_dir: str):
        pass
