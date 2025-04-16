from typing import Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pegasusio import MultimodalData
import pegasus as pg

from ddqc.filtering import perform_ddqc
from ddqc.plotting import filtering_facet_plot, boxplot_sorted, calculate_filtering_stats
from ddqc.utils import cluster_data, calculate_percent_ribo

#main changes: passing in a clustering, assuming intial QC done - because the point is not to have an end to end thing
#literally only need the mad calculation realistically
#question is do we want to generalize this to set of qc metrics that we want?
def ddqc_metrics(data: MultimodalData,
                 clustering_obs: str = "louvain",
                 default_threshold: float = 2.0, 
                 metrics_df: pd.Dataframe = None,
                 return_df_qc: bool = False, display_plots: bool = True) -> Union[None, pd.DataFrame]:
    """
    Parameters:
        data (MultimodalData): Pegasus object.
        clustering_obs (str): name of column in obs to use for clustering.
        default_threshold (float): parameter for the selected method (default is 2).
            Note that "outlier" method doesn't requre parameter and will ignore this option.
        metrics_df (dataframe): dataframe of metrics, index == metrics, columns store threshhold and bound information, 
            as well as which bounds to calculate
        return_df_qc (bool): whether to return a dataframe with the information about on what metric and what threshold
            the cell was removed for each removed cell. (default is False)
        display_plots (bool): whether to show plots that would show filtering statistics (default is True).
    Returns:
        (None, pandas.DataFrame). DataFrame with cluster labels and thresholds for each metric is returned
            if return_df_qc was True.
    """
    assert isinstance(data, MultimodalData)
    if 
    # initial qc
    #consider filling in percent calculations here
    passed_qc, df_qc, _ = perform_ddqc(data, clustering_obs, default_threshold, metrics_df)
    #TODO: below here will still need help.
    if display_plots:
        boxplot_sorted(df_qc, "n_genes", "cluster_labels", hline_x=np.log2(200), log=True)
        plt.show()
        boxplot_sorted(df_qc, "percent_mito", "cluster_labels", hline_x=10)
        plt.show()
        if ((threshold_counts == 0 or threshold_counts is None)
                and (threshold_genes == 0 or threshold_genes is None)
                and (threshold_mito == 0 or threshold_mito is None)
                and (threshold_ribo == 0 or threshold_ribo is None)):
            if method == "mad":
                fs = calculate_filtering_stats(data_copy, threshold, n_genes_lower_bound, percent_mito_upper_bound)
                filtering_facet_plot(fs, threshold, pct=False)
                plt.show()
                filtering_facet_plot(fs, threshold, pct=True)
                plt.show()

    # reverse_to_raw_matrix(data.current_data(), obs_copy, var_copy, uns_copy)
    data.obs["passed_qc"] = passed_qc
    df_qc["passed_qc"] = data.obs["passed_qc"]

    if return_df_qc:
        return df_qc
