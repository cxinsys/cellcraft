import os
import multiprocessing.pool as mpp

try:
    import torch
except (ModuleNotFoundError, ImportError) as err:
    pass

import numpy as np

def save_results(droot, rss, W, A, B, node_name):
    os.makedirs(droot, exist_ok=True)
    np.savetxt(os.path.join(droot, "RSS.txt"), [rss], delimiter="\t", fmt="%.14f")
    np.savetxt(os.path.join(droot, "W.txt"), W, delimiter="\t", fmt="%.14f")

    tmp_rm = np.concatenate([node_name[:, None], A.astype(str)], axis=1)
    extended_nn = np.concatenate((['TE'], node_name))
    tmp_rm = np.concatenate([extended_nn[None, :], tmp_rm])

    np.savetxt(os.path.join(droot, "A.txt"), tmp_rm, delimiter="\t", fmt="%s")
    np.savetxt(os.path.join(droot, "B.txt"), B, delimiter="\t", fmt="%.14f")


def calculate_lm_memory_usage(batch, exp_data_shape, new_b_shape, dtype=np.float32):
    # 데이터 타입 크기 (float64는 8바이트)
    dtype_size = np.dtype(dtype).itemsize

    # 크기 정보
    g, c = exp_data_shape  # exp_data 크기 (g: gene, c: cell)
    s, z = new_b_shape     # new_b 크기 (s: sampling batch, z: latent factor)

    # 입력 데이터 메모리 크기
    exp_data_memory = batch * c * dtype_size  # exp_data
    pseudotime_memory = c * dtype_size    # pseudotime
    new_b_memory = s * z * dtype_size         # new_b

    # 중간 계산 배열 메모리 크기
    Z_memory = s * z * c * dtype_size                     # Z
    XtX_memory = s * z * z * dtype_size                   # XtX
    Xty_memory = s * z * batch * dtype_size                   # Xty
    W_memory = s * batch * z * dtype_size                     # W
    WZ_memory = s * batch * c * dtype_size                    # WZ
    diffs_memory = s * batch * c * dtype_size                 # diffs

    fix_memory = (pseudotime_memory + new_b_memory + Z_memory + XtX_memory) / (1024 ** 2)
    other = (exp_data_memory + Xty_memory + W_memory + WZ_memory + diffs_memory + diffs_memory) / (1024 ** 2)

    return fix_memory, other


def check_gpu_computability(w_shape, new_b_shape, dtype=np.float32):
    # 데이터 타입 크기 (float64는 8바이트)
    dtype_size = np.dtype(dtype).itemsize

    # 크기 정보
    g, z = w_shape  # exp_data 크기 (g: gene, c: cell)
    z = new_b_shape[-1]     # new_b 크기 (z: latent factor)

    # 입력 데이터 메모리 크기
    b_matrix_memory = z * z * dtype_size  # exp_data
    w_memory = g * z * dtype_size    # pseudotime

    # 중간 계산 배열 메모리 크기
    inv_w_memory = z * g * dtype_size                     # Z

    total_memory = (b_matrix_memory + w_memory + inv_w_memory) / (1024 ** 2)

    gpu_mem = get_gpu_memory(gpu_index=0)

    return total_memory < gpu_mem


def get_gpu_memory(gpu_index=0):
    torch.cuda.set_device(gpu_index)  # GPU 선택
    total_memory = torch.cuda.get_device_properties(gpu_index).total_memory

    return 0.9 * (total_memory / (1024 ** 2))


def calculate_batchsize(batch,
                        exp_data_shape,
                        new_b_shape,
                        dtype=np.float64,
                        num_gpus=1,
                        num_ppd=1):
    fix, other = calculate_lm_memory_usage(batch=batch,
                                           exp_data_shape=exp_data_shape,
                                           new_b_shape=new_b_shape,
                                           dtype=dtype)
    gpu_mem = get_gpu_memory(gpu_index=0)

    free_mem = gpu_mem - (fix * num_ppd)

    if free_mem < 0:
        raise ValueError("The number of processors you want to use is too many for the batch size. ")


    remain_mem = free_mem - other

    dtype_size = np.dtype(dtype).itemsize

    one_batch_mem = exp_data_shape[1] * dtype_size / (1024 ** 2)
    outer_batch = np.ceil(remain_mem / (one_batch_mem * num_gpus)).astype(np.int32)

    if outer_batch < 0:
        raise ValueError("The number of processors you want to use is too many for the batch size. ")

    return outer_batch

def istarmap(self, func, iterable, chunksize=1):
    """starmap-version of imap
    """
    self._check_running()
    if chunksize < 1:
        raise ValueError(
            "Chunksize must be 1+, not {0:n}".format(
                chunksize))

    task_batches = mpp.Pool._get_tasks(func, iterable, chunksize)
    result = mpp.IMapIterator(self)
    self._taskqueue.put(
        (
            self._guarded_task_generation(result._job,
                                          mpp.starmapstar,
                                          task_batches),
            result._set_length
        ))
    return (item for chunk in result for item in chunk)


mpp.Pool.istarmap = istarmap