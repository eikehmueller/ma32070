program test_comm

  implicit none
  include "mpif.h"

  ! Minimal message size = 2**imin
  integer, parameter :: imin = 0
  ! Maximal message size = 2**imax
  integer, parameter :: imax = 24
  ! Number of sends/recvs
  integer, parameter :: nsend = 1000

  integer :: myid, numprocs, ierr
  ! Note these fancy dimension statements that allow a custom index range
  real(kind=8), dimension(imin:imax) :: t_elapsed
  integer, dimension(imin:imax) :: messagesize
  integer :: i

  ! Initialise MPI and work out number of processes and local rank
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr)

  ! Check that we are running on two processors
  if (numprocs /= 2) then
     if (myid == 0) then
        ! Often we want to print from only one process to avoid floods
        print *, "Program can only run on two processors. Aborting..."
     end if
     call MPI_Finalize(ierr)
     STOP
  end if

  if (myid == 0) then
     write(*,*) '*** Pingpong test ***'
  endif

  ! Loop over all message sizes
  do i = imin, imax
     messagesize(i) = 2**i     
     call pingpong(messagesize(i),nsend,t_elapsed(i))
     if (myid == 0) then
        write(*,'("Message size = ",I10," t = ",F12.4," micro seconds")') &
             messagesize(i), 1.E6*t_elapsed(i)
     end if
  end do

  ! Write results to disk in a format compatible with lsprogs
  if (myid == 0) then
     open(unit=2,file='lsdata.txt')
     write(2,*) imax-imin+1
     write(2,*) messagesize(imin:imax)
     write(2,*) t_elapsed(imin:imax)
     close(2)
  end if

    ! Finalise MPI
  call MPI_Finalize(ierr)

end program test_comm
