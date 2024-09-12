import os
import os.path as osp
import sys
import numpy as np
import fasttenet as fte

if __name__ == "__main__":
    exp_matrix = sys.argv[1]
    pseudotime = sys.argv[2]
    cell_select = sys.argv[3]
    result_matrix = sys.argv[4]
    backend = sys.argv[5]
    num_devices = int(sys.argv[6])
    batch_size = int(sys.argv[7])

    # Create worker
    worker = fte.FastTENET(dpath_exp_data=exp_matrix,  # Required
                           dpath_trj_data=pseudotime,  # Required
                           dpath_branch_data=cell_select,  # Required
                           spath_result_matrix=result_matrix,  # Optional
                           make_binary=True)  # Optional, default: False

    result_matrix = worker.run(backend=backend,
                               device_ids=num_devices,
                               procs_per_device=1,
                               batch_size=batch_size,  # k1 - 2080ti: 2**15, 3090: 2**16 / k3 - 2**14, 2**15
                               num_kernels=1,
                               method='shift_left',
                               kp=0.5)
