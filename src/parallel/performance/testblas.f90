program testblas

  implicit none

  ! compares times for lu1blas2 and lu1blas3 on random well-conditioned
  ! test matrices
 
  integer :: n
  real (kind=8),  dimension(:,:), allocatable :: A, Aorig
  real (kind=8),  dimension(:), allocatable :: b,x_e

  real :: t1, t2, t_blas2, t_blas3, t_flop, t_mem
  integer :: bsize

  n = 4096
  print *, 'n = ', n

  ! allocate storage
  allocate(A(n,n), Aorig(n,n), x_e(n), b(n))

  ! create a random well-conditioned symmetric positive definite matrix
  call rsymm(A,n)
  Aorig = A
 
  !  random solution vector
  call random_number(x_e) 
  ! assemble RHS: so that Ax = b has the exact solution x_e
  b  = matmul(A,x_e)

  !!!!! LU1 using BLAS up to Level 2 subroutines  
  print *, 'lublas2 started ...'
  call cpu_time(t1)
  call lublas2(A,n,n)
  call cpu_time(t2)
  t_blas2 = t2-t1
  print *, 'lublas2 took ', t_blas2, 'seconds'
  call check_lu(n,A,b,x_e)
  A = Aorig
    
  !!!!! LU1 using BLAS Level 3
  ! use fixed block size up to 128 for the moment
  bsize = min(n,128)
  print *, 'lublas3 started ...'
  call cpu_time(t1)
  call lublas3(A,n,bsize)
  call cpu_time(t2)
  t_blas3 = t2-t1
  print *, 'lublas3 took ', t_blas3, 'seconds'
  call check_lu(n,A,b,x_e)

  deallocate(A,Aorig,b,x_e)

  t_flop = t_blas3 / (2./3.*(1.*n)**3)
  t_mem = t_blas2 / (2./3.*(1.*n)**3) - t_flop

  write(*,'(A16,E10.4,A2)') "t_flop = ",t_flop,"s"
  write(*,'(A16,E10.4,A2)') "t_mem = ",t_mem,"s"

end program testblas
