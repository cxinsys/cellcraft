import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
fastscode_root = os.path.abspath(os.path.join(current_dir, 'fastscode'))
sys.path.append(current_dir)

import os.path as osp

import numpy as np

from fastscode.fastscode import FastSCODE
from fastscode.inference.inference import NetWeaver



if __name__ == "__main__":
    dpath_exp_data = sys.argv[1]  # expression data
    dpath_trj_data = sys.argv[2]  # pseudotime data
    dpath_branch_data = sys.argv[3]  # cell select data

    num_z = int(sys.argv[4])  # number of z
    max_iter = int(sys.argv[5])  # number of iteration

    backend = sys.argv[6]  # backend, 'gpu', 'torch', etc..
    num_devices = int(sys.argv[7])  # number of devices
    sb = int(sys.argv[8])  # sampling batch size
    chunk_size = int(sys.argv[9])  # chunk size

    exp_data = np.loadtxt(dpath_exp_data, delimiter=",", dtype=str)  # cell x gene
    node_name = exp_data[0, 1:]
    exp_data = exp_data[1:, 1:].astype(np.float64).T  # gene x cell
    pseudotime = np.loadtxt(dpath_trj_data, delimiter="\t")
    branch = np.loadtxt(dpath_branch_data, delimiter="\t")

    pseudotime = pseudotime[branch == 1]
    exp_data = exp_data[:, branch == 1]

    repeats = np.arange(int(sys.argv[10]), dtype=np.int32)  # num repeat

    links = int(sys.argv[11])
    trim_threshold = float(sys.argv[12])

    As = []
    for r in repeats:
        worker = FastSCODE(exp_data=exp_data,
                              pseudotime=pseudotime,
                              node_name=node_name,
                              droot=None,
                              num_tf=None,
                              num_cell=None,
                              num_z=num_z,
                              max_iter=max_iter,
                              dtype=np.float64)

        rss, W, A, B = worker.run(backend=backend,
                                  device_ids=num_devices,
                                  sampling_batch=sb,
                                  chunk_size=chunk_size)

        As.append(A)

    mean_A = np.mean(As, axis=0)
    # RSS, W, A and B will save under the droot directory

    weaver = NetWeaver(result_matrix=mean_A,
                                   gene_names=node_name,
                                   links=links,
                                   is_trimming=True,
                                   trim_threshold=trim_threshold,
                                   dtype=np.float64)

    grns = weaver.run(backend=backend, device_ids=num_devices, batch_size=100)

    dpath_grn = osp.abspath(osp.dirname(dpath_exp_data))

    grn, trimmed_grn = grns
    outdegrees = weaver.count_outdegree(grn)
    trimmed_ods = weaver.count_outdegree(trimmed_grn)

    fpath_save_grn = sys.argv[13]  # grn file name
    np.savetxt(fpath_save_grn, grn, delimiter='\t', fmt="%s")
    print('save grn in ', fpath_save_grn)

    fpath_save_odg = sys.argv[14]  # outdegree file name
    np.savetxt(fpath_save_odg, outdegrees, fmt="%s")
    print('save grn outdegrees in ', fpath_save_odg)

    fpath_save_trimmed_grn = sys.argv[15]  # trimmed grn file name
    np.savetxt(fpath_save_trimmed_grn, trimmed_grn, delimiter='\t', fmt="%s")
    print('save trimmed grn in ', fpath_save_trimmed_grn)

    fpath_save_trimmed_odg = sys.argv[16]  # trimmed outdegree file name
    np.savetxt(fpath_save_trimmed_odg, trimmed_ods, fmt="%s")
    print('save trimmed grn outdegrees in ', fpath_save_trimmed_odg)
