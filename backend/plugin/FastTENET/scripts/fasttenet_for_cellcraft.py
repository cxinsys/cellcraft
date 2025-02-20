import os
import sys

import os.path as osp

import numpy as np

# from fasttenet.fasttenet import FastTENET
# from fasttenet.utils import *
# from fasttenet.inference.inference import NetWeaver
import fasttenet as fte

if __name__ == "__main__":

    dpath_exp_data = sys.argv[1]  # expression data
    dpath_trj_data = sys.argv[2]  # pseudotime data
    dpath_branch_data = sys.argv[3]  # cell select data
    dpath_tf_data = sys.argv[4]  # tf data

    fpath_save = sys.argv[5]
    fpath_timmed = sys.argv[6]
    fpath_out = sys.argv[7]
    fpath_timmed_out = sys.argv[8]

    spath_result_matrix = None

    backend = sys.argv[9]
    num_devices = int(sys.argv[10])
    batch_size = int(sys.argv[11])

    fdr = float(sys.argv[12])
    links = int(sys.argv[13])
    trim_threshold = float(sys.argv[14])

    node_name, exp_data = fte.load_exp_data(dpath_exp_data, make_binary=True)
    trajectory = fte.load_time_data(dpath_trj_data, dtype=np.float32)
    branch = fte.load_time_data(dpath_branch_data, dtype=np.int32)

    tf = None
    if dpath_tf_data != "None":
        tf = np.loadtxt(dpath_tf_data, dtype=str)
    # if spath_result_matrix == "None":
    #     spath_result_matrix = None

    aligned_data = fte.align_data(data=exp_data, trj=trajectory, branch=branch)

    # Create worker
    # expression data, trajectory data, branch data path is required
    # tf data path is optional
    # save path is optional
    worker = fte.FastTENET(aligned_data=aligned_data, # Required
                           node_name=node_name, # Required
                           tfs=tf, # Required
                           spath_result_matrix=spath_result_matrix, # Optional
                           make_binary=True) # Optional, default: False

    result_matrix = worker.run(backend=backend, device_ids=num_devices, procs_per_device=1, batch_size=batch_size,
                               num_kernels=1, binning_method='FSBW-L', kp=0.5)

    weaver = fte.NetWeaver(result_matrix=result_matrix,
                                   gene_names=node_name,
                                   tfs=tf,
                                   fdr=fdr,
                                   links=links,
                                   is_trimming=True,
                                   trim_threshold=trim_threshold,
                                   dtype=np.float32
                                   )

    grn, trimmed_grn = weaver.run(backend=backend,
                                  device_ids=num_devices,
                                  batch_size=0)

    outdegrees = weaver.count_outdegree(grn)
    trimmed_ods = weaver.count_outdegree(trimmed_grn)

    dpath_grn = osp.abspath(osp.dirname(dpath_exp_data))
    # if links != 0:
    #     fpath_save = osp.join(dpath_grn, f"result_matrix.links" + str(links) + ".sif")
    # else:
    #     fpath_save = osp.join(dpath_grn, f"result_matrix.fdr" + str(fdr) + ".sif")

    np.savetxt(fpath_save, grn, delimiter='\t', fmt="%s")
    print('save grn in ', fpath_save)

    np.savetxt(fpath_timmed, trimmed_grn, delimiter='\t', fmt="%s")
    print('save trimmed grn in ', fpath_timmed)

    np.savetxt(fpath_out, outdegrees, fmt="%s")
    print('save grn outdegrees in ', fpath_out)

    np.savetxt(fpath_timmed_out, trimmed_ods, fmt="%s")
    print('save trimmed grn outdegrees in ', fpath_timmed_out)