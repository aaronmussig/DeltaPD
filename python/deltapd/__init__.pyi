from enum import Enum

class PyDeltaPD:

    def __init__(self, qry_dm: PyDistMatrix, ref_dm: PyDistMatrix):
        pass

    def run(self):
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
    def __init__(self, knn: int, cpus: int, model: PyLinearModelType, error: PyLinearModelError, corr: PyLinearModelCorr):
        pass


