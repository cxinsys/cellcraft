import multiprocessing

from tqdm import tqdm
import numpy as np
import torch

from mate.array import get_array_module
from mate.utils import get_device_list


class Simulator(object):
    def __init__(self,
                 init=None,
                 result_matrix=None,
                 num_times=100,
                 fpath_out=None):

        self.init = init
        self.rm = result_matrix
        self.tmp_t = np.linspace(0, 1, num_times + 1)
        self.fpath_out = fpath_out

    def run(self,
            backend=None,
            device_ids=None,
            batch_size=0):

        if not backend:
            backend = 'cpu'

        if not device_ids:
            if backend == 'cpu':
                device_ids = [0]
            else:
                device_ids = get_device_list()

        if type(device_ids) is int:
            list_device_ids = [x for x in range(device_ids)]
            device_ids = list_device_ids

        if not batch_size:
            batch_size = 10

        print("Computes the eigenvalue...")
        l, U = np.linalg.eig(self.rm)
        print("Done")

        print("[DEVICE: {}, Num. Processor: {}, Batch Size: {}]".format(backend, len(device_ids), batch_size))

        multiprocessing.set_start_method('spawn', force=True)

        with multiprocessing.Pool(processes=len(device_ids)) as pool:
            list_backend = []
            list_id = []
            list_init = []
            list_l = []
            list_U = []
            list_tmp_t = []
            list_batch = []

            outer_batch = np.ceil(len(self.tmp_t) / (len(device_ids))).astype(np.int32)
            for i, start in enumerate(range(0, len(self.tmp_t), outer_batch)):
                end = start + outer_batch

                list_backend.append(backend + ":" + str(device_ids[i % len(device_ids)]))
                list_id.append(i)
                list_init.append(self.init)
                list_l.append(l)
                list_U.append(U)
                list_tmp_t.append(self.tmp_t[start:end])
                list_batch.append(batch_size)

            inputs = zip(list_backend, list_id, list_init, list_l, list_U, list_tmp_t, list_batch)

            list_mean = []

            for batch_result in pool.istarmap(self.compute, inputs):
                list_mean.extend(batch_result)

        result = np.concatenate(list_mean).T
        result = np.concatenate((self.tmp_t[None, :], result), axis=0)

        if self.fpath_out is not None:
            np.savetxt(self.fpath_out, result, delimiter="\t", fmt="%.14f")

        return result

    def compute(self,
                backend=None,
                id=0,
                init=None,
                l=None,
                U=None,
                tmp_t=None,
                batch_size=None,
                ):

        am = get_array_module(backend)

        init = am.array(init, dtype=str(init.dtype))
        tmp_t = am.array(tmp_t, dtype=str(tmp_t.dtype))
        l = am.array(l, dtype=str(l.dtype))
        U = am.array(U, dtype=str(U.dtype))

        invU = am.inv(U)

        list_mean = []
        for i , start in enumerate(tqdm(range(0, len(tmp_t), batch_size), desc=f"Process {id}", position=id, leave=True)):
            end = start + batch_size

            exp = am.exp(l[None, :] * tmp_t[start:end, None])
            diag = am.zeros((exp.shape[0], exp.shape[1], exp.shape[1]), dtype=str(exp.dtype))
            i = am.arange(exp.shape[1])
            diag[:, i, i] = exp

            at1 = am.dot(U, am.transpose(diag, axes=(1, 2, 0)))
            at2 = am.dot(am.transpose(at1, axes=(2, 0, 1)), invU)
            eAt = am.real(at2)
            mean = am.dot(eAt, init)

            list_mean.append(am.asnumpy(mean))

        return list_mean
