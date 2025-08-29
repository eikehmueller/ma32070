"""Driver for parallel matrix-vector multiplication

Evaluate performance of matrix-matrix multiplication in matmul.py for
different problem sizes and setups.

Run as

    python driver.py --help

to list command line options.
"""

import argparse
from mpi4py import MPI
import numpy as np
from matmul import parallel_matmul

parser = argparse.ArgumentParser()

parser.add_argument(
    "--n",
    action="store",
    default=128,
    required=True,
    type=int,
    help="Number of rows of matrix A",
)

parser.add_argument(
    "--m",
    action="store",
    required=False,
    type=int,
    help="Number of columns of matrix A",
)

parser.add_argument(
    "--r",
    action="store",
    required=False,
    type=int,
    help="Number of columns of matrix B",
)

parser.add_argument(
    "--overlap",
    action="store_true",
    help="Overlap calculation and communication?",
)

parser.add_argument(
    "--nocomm",
    action="store_true",
    help="No parallel communication (leads to incorrect results)?",
)

parser.add_argument(
    "--nocheck",
    action="store_true",
    help="Do not check correctness of results?",
)

args, _ = parser.parse_known_args()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nproc = comm.Get_size()

n = args.n
m = args.m if args.m else args.n
r = args.r if args.r else args.n

assert n % nproc == 0, "n must be a multiple of the number of processors"
assert m % nproc == 0, "m must be a multiple of the number of processors"

n_loc = n // nproc
m_loc = m // nproc

if rank == 0:
    print(f"Running on {nproc} processors.")
    print(f" n = {n:8d}")
    print(f" m = {m:8d}")
    print(f" r = {r:8d}")
    print(f" overlap = {args.overlap}")
    print(f" nocomm  = {args.nocomm}")
    print(f" nocheck = {args.nocheck}")
    print()
rng = np.random.default_rng(seed=411587)
if not args.nocheck:
    # Global matrices
    A = rng.normal(size=(n, m))
    B = rng.normal(size=(m, r))
    # Result for testing
    C_true = A @ B
t_elapsed = 0
niter = 0
t_target = 10.0
while t_elapsed < t_target:
    if args.nocheck:
        A_loc = rng.normal(size=(n_loc, m))
        B_loc = rng.normal(size=(m_loc, r))
        C_loc = np.zeros(shape=(n_loc, r))
    else:
        A_loc = np.array(A[rank * n_loc : (rank + 1) * n_loc, :])
        B_loc = np.array(B[rank * m_loc : (rank + 1) * m_loc, :])
        C_loc = np.zeros(shape=(n_loc, r))
    comm.Barrier()
    t_start = MPI.Wtime()
    parallel_matmul(A_loc, B_loc, C_loc, args.overlap, args.nocomm)
    comm.Barrier()
    t_finish = MPI.Wtime()
    delta_t = comm.allreduce(t_finish - t_start, MPI.MAX)
    t_elapsed += delta_t
    niter += 1

t_elapsed /= niter

if not args.nocheck:
    C = np.zeros(shape=(n, r))
    comm.Allgather(C_loc, C)
    nrm = np.linalg.norm(C - C_true) / np.linalg.norm(C_true)

if rank == 0:
    if not args.nocheck:
        print(f"error norm = {nrm:8.4e}")
    print(f"elapsed time = {t_elapsed:8.2e} s")
