from enum import Enum

from deltapd import PyLinearModelError, PyLinearModelCorr


class CorrelationFn(str, Enum):
    R2 = "R2"
    Pearson = "Pearson"

    def to_rs(self):
        if self is CorrelationFn.R2:
            return PyLinearModelCorr.R2
        elif self is CorrelationFn.Pearson:
            return PyLinearModelCorr.Pearson
        else:
            raise ValueError(f"Unknown correlation function: {self}")


class ErrorFn(str, Enum):
    MSE = "MSE"
    NormMSE = "NormMSE"
    RMSE = "RMSE"

    def to_rs(self):
        if self is ErrorFn.MSE:
            return PyLinearModelError.MSE
        elif self is ErrorFn.NormMSE:
            return PyLinearModelError.NormMSE
        elif self is ErrorFn.RMSE:
            return PyLinearModelError.RMSE
        else:
            raise ValueError(f"Unknown error function: {self}")
