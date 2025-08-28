subroutine rsymm(A,n)

  implicit none

  ! assembles a random well-conditioned n times n symmetric matrix A
  ! It is just a small symmetric random perturbation of the identity
  ! and so should be well-conditioned

  integer, intent(in) :: n
  real (kind=8), dimension(n,n), intent(out) :: A
  integer :: i
  integer,dimension(129) :: seed = 0

  call random_seed(put=seed)
  call random_number(A)

  A = 1d-2*(A*transpose(A))    ! small symmetric matrix

  do i = 1,n
    A(i,i) = A(i,i) + 1.0_8        
  end do

end subroutine rsymm
