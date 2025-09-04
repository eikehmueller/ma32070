import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

data = np.zeros(4, dtype=int)
data[:] = rank + 1

print(f"BEFORE: rank {rank} has data", data)

if rank == 0:
    comm.Send(data, 1)
else:
    comm.Recv(data, source=0)

print(f"AFTER: rank {rank} has data", data)
