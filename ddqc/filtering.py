from typing import Union

import numpy as np
import pandas as pd
from pegasusio import MultimodalData

from ddqc.utils import mad


def metric_filter(data: MultimodalData, param: float, metric_name: str, metric_info: pd.Series,
                  df_qc: pd.DataFrame = None) -> np.ndarray:
    """Function that performs filtering result computation for the specified metric."""
    param   =  m_info["threshold"] if m_info["threshold"] is not None else default_threshold
    lower_bound = m_info["lower_bound"]
    upper_bound = m_info["upper_bound"]
    qc_pass = np.zeros(data.shape[0], dtype=bool)  # T/F array to tell whether the cell is filtered
    if df_qc is not None:
        df_qc[f"{metric_name}_lower_co"] = None
        df_qc[f"{metric_name}_upper_co"] = None

    for cl in data.obs["cluster_labels"].cat.categories:  # iterate through all clusters
        idx = data.obs["cluster_labels"] == cl
        values = data.obs.loc[idx, metric_name]

        # calculate MAD cutoffs, which are median ± param * MAD
        median_v = np.median(values)
        mad_v = mad(values)
        lower_co = median_v - param * mad_v
        upper_co = median_v + param * mad_v

        lower_co = min(lower_co, lower_bound) if lower_bound is not None else lower_co
        upper_co = max(upper_co, upper_bound) if upper_bound is not None else upper_co

        qc_pass_cl = np.ones(values.size, dtype=bool)
        if df_qc is not None:
            df_qc.loc[idx, f"{metric_name}"] = values
        if m_info["do_lower"]:
            qc_pass_cl &= (values >= lower_co)
            if df_qc is not None:
                df_qc.loc[idx, f"{metric_name}_lower_co"] = lower_co
        if m_info["do_upper"]:
            qc_pass_cl &= (values <= upper_co)
            if df_qc is not None:
                df_qc.loc[idx, f"{metric_name}_upper_co"] = upper_co
        if df_qc is not None:
            df_qc.loc[idx, f"{metric_name}_passed_qc"] = qc_pass_cl
        qc_pass[idx] = qc_pass_cl
    return qc_pass


def perform_ddqc(data: MultimodalData, clustering_obs: str,
                 default_threshold: float, 
                 metrics_df: pd.Dataframe,
                 filtering_stats: Union[pd.DataFrame, None] = None):
    """Function that computes ddqc for all requested metrics with specified parameters."""
    df_qc = pd.DataFrame({clustering_obs: data.obs[clustering_obs].values}, index=data.obs_names)
    passed_qc = np.ones(data.shape[0], dtype=bool)
    #question is how we want to deal with, do we want a data frame and what passes, could see reason, but do we need at cell level?
    #going to end up being a space time tradeoff... maybe we do accept this for sake of sanity

    for me, m_info in metrics_df.iterrows():
      passed_qc &= metric_filter(data, default_threshold, metric, m_info, df_qc)

    if filtering_stats is not None:
        cl_filtering_stats = []
        for cl in data.obs[clustering_obs].cat.categories:
            idx = data.obs[clustering_obs] == cl
            unique, counts = np.unique(passed_qc[idx], return_counts=True)
            n_filtered_cells = dict(zip(unique, counts)).get(False, 0)
            n_filtered_cells_pct = n_filtered_cells / len(passed_qc[idx]) * 100
            row = {"threshold": threshold, "cluster": cl, "filtered_cells": n_filtered_cells,
                   "filtered_cells_pct": n_filtered_cells_pct}
            cl_filtering_stats.append(row)
        filtering_stats = pd.concat([filtering_stats, pd.DataFrame.from_records(cl_filtering_stats)])
    return passed_qc, df_qc, filtering_stats
