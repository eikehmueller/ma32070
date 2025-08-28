subroutine lublas2(A,n,lda)
  implicit none
  ! Gaussian Elimination without pivoting
  ! using blas level 2 for a top square n x n submatrix
  ! of a lda x n matrix A
  
  integer, intent(in) :: n, lda
  real (kind=8), intent(inout), dimension(lda,n) :: A
  integer :: i

  do i = 1,n-1
     ! scaling x = alpha*x
     call dscal(n-i   ,   1.0_8/A(i,i),  A(i+1,i),1          )
     !          length,   alpha       ,  x   &    stride of x
     
     ! Rank 1 Update A = A + alpha*v*w^T
     call dger(n-i,n-i, -1.0_8, A(i+1,i),1, A(i,i+1),lda, A(i+1,i+1),lda)
     !          n , m ,  alpha, v & stride, w^T & stride, A & lda
  end do

end subroutine lublas2
