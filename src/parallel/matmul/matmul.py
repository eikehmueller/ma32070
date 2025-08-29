"""Parallel matrix-vector multiplication

Computes

    C = A @ B

in parallel, where the shapes of the matrices are:

    * A: (n,m)
    * B: (m,r)
    * C: (n,r)

n and m have to be multiples of the number of processors.

Supports overlapping of computation and communications.
"""

from mpi4py import MPI
import numpy as np


def parallel_matmul(A, B, C, overlap=False, nocomm=False):
    """Parallel matrix-vector product

    Compute local part of matrix C = A @ B, where each processor only stores its
    local part of the matrices A and B.

    :arg A: process-local part of matrix A
    :arg B: process-local part of matrix B
    :arg C: process-local part of resulting matrix C
    :arg overlap: overlap computation and communication?

    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nproc = comm.Get_size()

    m_loc, _ = B.shape

    B_recv = np.empty_like(B)
    req_send = None
    req_recv = None
    for q in range(nproc):
        if not nocomm and overlap and q < nproc - 1:
            req_send = comm.Isend(B, dest=(rank - 1) % nproc)
            req_recv = comm.Irecv(B_recv)
        C[:, :] += (
            A[:, ((q + rank) % nproc) * m_loc : (((q + rank) % nproc + 1)) * m_loc]
            @ B[:, :]
        )

        if q < nproc - 1 and not nocomm:
            if overlap:
                req_send.wait()
                req_recv.wait()
                B[:, :] = B_recv[:, :]
            else:
                comm.Sendrecv_replace(B, (rank - 1) % nproc)
    return C
