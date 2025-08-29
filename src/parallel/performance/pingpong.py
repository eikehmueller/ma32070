import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD

rank = comm.Get_rank()
nproc = comm.Get_size()

assert nproc == 2, "Can only be run on exactly two ranks"

for n in 2 ** np.arange(24, dtype=int):
    recv_buf = np.ones(n)
    send_buf = np.ones(n)
    t_elapsed = 0
    niter = 0
    comm.barrier()
    while t_elapsed < 0.1:
        t_start = MPI.Wtime()
        if rank == 0:
            comm.Send(send_buf, 1)
            comm.Recv(recv_buf, source=1)
        else:
            comm.Recv(recv_buf, source=0)
            comm.Send(send_buf, 0)
        t_finish = MPI.Wtime()
        t_elapsed += comm.allreduce(t_finish - t_start, MPI.MAX)
        niter += 1
    comm.barrier()
    t_elapsed /= 2 * niter
    if rank == 0:
        print(f"Python message size = {n:10d} t = {1E6*t_elapsed:12.4} micro seconds")
