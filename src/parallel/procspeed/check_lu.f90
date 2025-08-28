subroutine check_lu(n,A,b,x_e)
  implicit none
  integer, intent(in) :: n
  real(kind=8), dimension(n,n), intent(in) :: A
  real(kind=8), dimension(n), intent(in) :: b,x_e

  ! local vector
  real(kind=8), dimension(n) :: x

  
  ! Back substitution
  call bsblas2(n,A,b,x)
  
  ! check that the solver has worked
  print*, 'error in solve: ', norm2(x-x_e)
  print*, ''

end subroutine check_lu
