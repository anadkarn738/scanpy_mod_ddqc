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
                 method: str = "mad", threshold: float = 2.0, threshold_counts: Union[int, None] = 0,
                 threshold_genes: Union[int, None] = 0, threshold_mito: Union[float, None] = 0,
                 threshold_ribo: Union[float, None] = 0,
                 mito_prefix: str = "MT-", ribo_prefix: str = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",
                 n_genes_lower_bound: int = 200, percent_mito_upper_bound: float = 10.0, random_state: int = 29,
                 return_df_qc: bool = False, display_plots: bool = True) -> Union[None, pd.DataFrame]:
    """
    Parameters:
        data (MultimodalData): Pegasus object.
        clustering_obs (str): name of column in obs to use for clustering.
        method (str): statistic on which the threshold would be calculated. Supported options are "mad" and "outlier"
            (default is "mad").
        threshold (float): parameter for the selected method (default is 2).
            Note that "outlier" method doesn't requre parameter and will ignore this option.
        threshold_counts (int, None): setting for applying ddqc based on number of counts. (Default is 0)
            - If set to 0, will perform ddqc on number of counts using the "threshold" parameter provided earlier.
            - If set to a number other than 0, will  overwrite "threshold" parameter for number of counts.
            - If set to None, won't perform ddqc on number of counts.
        threshold_genes (int, None): Same as above, but for number of genes.
        threshold_mito (float, None): Same as above, but for percent of mitochondrial transcripts.
        threshold_ribo (float, None): Same as above, but for percent of ribosomal transcripts.
        mito_prefix (str): gene prefix used to calculate percent_mito in a cell (default is "MT-").
        ribo_prefix (str): gene regular expression used to calculate percent_ribo in a cell
            (default is "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA").
        n_genes_lower_bound (int): bound for lower n_genes cluster-level threshold (default is 200).
        percent_mito_upper_bound (float): bound for upper percent_mito cluster-level threshold (default is 10).
        random_state (int): random seed for clustering results reproducibility (default is 29)
        return_df_qc (bool): whether to return a dataframe with the information about on what metric and what threshold
            the cell was removed for each removed cell. (default is False)
        display_plots (bool): whether to show plots that would show filtering statistics (default is True).
    Returns:
        (None, pandas.DataFrame). DataFrame with cluster labels and thresholds for each metric is returned
            if return_df_qc was True.
    """
    assert isinstance(data, MultimodalData)

    # initial qc
    #consider filling in percent calculations here
    calculate_percent_ribo(data, ribo_prefix)  # calculate percent ribo

    passed_qc, df_qc, _ = perform_ddqc(data_copy, method, threshold,
                                       threshold_counts, threshold_genes, threshold_mito, threshold_ribo,
                                       n_genes_lower_bound, percent_mito_upper_bound)

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
