subroutine bsblas2(n,A,b,x)
  implicit none  

  ! Back substitution to solve A*x = b
  ! assuming that A contains its LU factors 
  ! using level 2 blas

  integer, intent(in) :: n
  real (kind=8), intent(in), dimension(n,n) :: A
  real (kind=8), intent(in), dimension(n) :: b
  real (kind=8), intent(out), dimension(n) :: x

  ! avoid overwriting b:
  x = b

  ! first solve L*y = b  Writes the answer into x
  ! L is Lower Unit triangular, stored Not transposed
  call dtrsv('Lower','Not transp','Unit',n,A,n,x,1)  

  ! Now solve U*x = y   again writes the answer into x
  call dtrsv('Upper','Not transp','Non-unit',n,A,n,x,1)

end subroutine bsblas2
