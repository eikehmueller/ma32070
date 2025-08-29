subroutine pingpong(messagesize, nsend, t_elapsed)

  implicit none
  include "mpif.h"
  
  ! Size of message to exchange
  integer, intent(in) :: messagesize 
  ! Number of sends/recvs
  integer, intent(in) :: nsend
  ! Elapsed time for one send/recv
  real(kind=8), intent(out) :: t_elapsed
  ! Send and receive buffers
  real(kind=8), dimension(messagesize,nsend) :: buffer_send, buffer_recv
  integer :: i, ierr, myid
  integer :: istat(MPI_STATUS_SIZE)
  real(kind=8) :: t_start, t_finish

  call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  call random_number(buffer_send)

  ! Start the clock
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  t_start = MPI_Wtime()

  ! Exchange message nsend times
  do i = 1, nsend
     if (myid == 0) then
        ! Process 0 sends first, then receives
        call MPI_Send(buffer_send(:,i), messagesize, &
             MPI_DOUBLE_PRECISION, 1, 0, MPI_COMM_WORLD, ierr)
        call MPI_Recv(buffer_recv(:,i), messagesize, &
             MPI_DOUBLE_PRECISION, 1, 1, MPI_COMM_WORLD, istat, ierr)
     else 
        ! Process 1 receives first, then sends
        call MPI_Recv(buffer_recv(:,i), messagesize, &
             MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, istat, ierr)
        call MPI_Send(buffer_send(:,i), messagesize, &
             MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, ierr)
     end if
  end do

  ! Stop the clock
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  t_finish = MPI_Wtime()

  t_elapsed = (t_finish-t_start)/(2.0_8*nsend)

end subroutine pingpong
