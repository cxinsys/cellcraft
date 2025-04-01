import os
import sys
import os.path as osp

import numpy as np
import pandas as pd
from arboreto.algo import genie3
from arboreto.utils import load_tf_names
from distributed import Client, LocalCluster
import scipy.stats
import statsmodels.sandbox.stats.multicomp
import networkx as nx

if __name__ == "__main__":
    dpath_exp_data = sys.argv[1]  # expression file path
    dpath_trj_data = sys.argv[2]  # pseudotime file path
    dpath_branch_data = sys.argv[3]  # cell select file path

    dpath_tf_data = None if sys.argv[4] == "None" else sys.argv[4] # optional, transcript factor file path

    # If both parameters are 0 (default), it outputs the full GRN.
    fdr = float(sys.argv[5])  # OPTIONAL, default 0, min: 0, max: 1
    links = int(sys.argv[6])  # OPTIONAL, default: 0, min: 0

    fpath_save_grn = sys.argv[7]  # output file name of grn
    fpath_save_odg = sys.argv[8]  # output file name of outdegrees

    expdata = pd.read_csv(dpath_exp_data, header=0, index_col=0)
    pseudotime = np.loadtxt(dpath_trj_data, delimiter="\t")
    branch = np.loadtxt(dpath_branch_data, delimiter="\t")

    tfs = 'all' if dpath_tf_data is None else load_tf_names(dpath_tf_data)

    pseudotime = pseudotime[branch == 1]
    expdata = expdata[branch == 1]  # cell x gene

    inds_sorted_trj = np.argsort(pseudotime)
    aligned_data = expdata.iloc[inds_sorted_trj]

    client = Client(processes=False)

    network = genie3(aligned_data.to_numpy(), client_or_address=client, gene_names=aligned_data.columns, tf_names=tfs)
    # source, target, value,

    # cols = list(network.columns)
    # cols[1], cols[2] = cols[2], cols[1]

    # network = network[cols]

    te = network.iloc[:, 2].to_numpy()
    source = network.iloc[:, 0].to_numpy()
    target = network.iloc[:, 1].to_numpy()

    if fdr != 0:
        te_zscore = (te - np.mean(te)) / np.std(te)
        te_pval = 1 - scipy.stats.norm.cdf(te_zscore)
        te_fdr = statsmodels.sandbox.stats.multicomp.multipletests(te_pval, alpha=0.05, method='fdr_bh')

        inds_cutoff = te_fdr[1] < fdr

        source = source[inds_cutoff]
        target = target[inds_cutoff]
        te = te[inds_cutoff]
    elif links != 0:
        sorted_inds = np.argsort(te)
        te = te[sorted_inds][::-1][:links]
        source = source[sorted_inds][::-1][:links]
        target = target[sorted_inds][::-1][:links]

    te_grn = np.stack((source, te, target), axis=1)

    dg = nx.from_edgelist(te_grn[:, [0, 2]], create_using=nx.DiGraph)
    outdegrees = sorted(dg.out_degree, key=lambda x: x[1], reverse=True)

    np.savetxt(fpath_save_grn, te_grn, delimiter='\t', fmt="%s")
    np.savetxt(fpath_save_odg, outdegrees, fmt="%s")