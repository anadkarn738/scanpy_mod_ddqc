import numpy as np
import pandas as pd
import pegasus as pg
from pegasusio import UnimodalData, MultimodalData


def mad(x: np.ndarray, constant: float = 1.4826) -> float:
    """Function that computes adjusted MAD for numpy array"""
    return constant * np.median(np.absolute(x - np.median(x)))

def reverse_to_raw_matrix(unidata: UnimodalData, obs_copy: pd.DataFrame, var_copy: pd.DataFrame, uns_copy: dict):
    """Function that reverses a Pegasus object to a raw matrix and removes all additional information"""
    unidata.obs = obs_copy
    unidata.var = var_copy
    unidata.matrices["X"] = unidata.matrices.pop("counts")
    unidata.obsm.clear()
    unidata.varm.clear()
    unidata.uns = uns_copy
