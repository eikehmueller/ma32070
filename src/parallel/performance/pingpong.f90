subroutine pingpong(messagesize, t_elapsed)

  implicit none
  include "mpif.h"
  
  ! Size of message to exchange
  integer, intent(in) :: messagesize 
  ! Elapsed time for one send/recv
  real(kind=8), intent(out) :: t_elapsed
  ! Send and receive buffers
  integer, parameter :: nmaxsend = 10000
  real(kind=8), dimension(messagesize,0:nmaxsend) :: buffer_send, buffer_recv
  integer :: ierr, myid
  integer :: istat(MPI_STATUS_SIZE)
  integer :: nsend
  real(kind=8) :: t_start, t_finish, t_delta, t_delta_global

  call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  call random_number(buffer_send)

  ! Start the clock
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  nsend = 0
  t_elapsed = 0
  ! Exchange message nsend times
  do while ((t_elapsed < 0.1) .and. (nsend < nmaxsend))
     t_start = MPI_Wtime()
     if (myid == 0) then
        ! Process 0 sends first, then receives
        call MPI_Send(buffer_send(:,nsend), messagesize, &
             MPI_DOUBLE_PRECISION, 1, 0, MPI_COMM_WORLD, ierr)
        call MPI_Recv(buffer_recv(:,nsend), messagesize, &
             MPI_DOUBLE_PRECISION, 1, 1, MPI_COMM_WORLD, istat, ierr)
     else 
        ! Process 1 receives first, then sends
        call MPI_Recv(buffer_recv(:,nsend), messagesize, &
             MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, istat, ierr)
        call MPI_Send(buffer_send(:,nsend), messagesize, &
             MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, ierr)
     end if
     t_finish = MPI_Wtime()
     t_delta = t_finish - t_start
     call MPI_Allreduce(t_delta,t_delta_global,1, &
          MPI_DOUBLE_PRECISION, MPI_MAX, &
          MPI_COMM_WORLD,ierr)
     t_elapsed = t_elapsed + t_delta_global
     nsend = nsend + 1
  end do

  ! Stop the clock
  call MPI_Barrier(MPI_COMM_WORLD,ierr)  

  t_elapsed = t_elapsed/(2.0_8*nsend)

end subroutine pingpong
