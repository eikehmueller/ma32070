subroutine lublas3(A,n,bsize)

  implicit none

  ! Gaussian Elimination without pivoting
  ! using blocking and blas level 3 for a square n X n matrix
  ! Assumes the blocksize divides n 
  
  integer, intent(in) :: n, bsize
  real (kind=8), intent(inout), dimension(n,n) :: A
  integer :: nb                        ! nb = number of blocks
  integer ::  i, j

  nb = n/bsize ! this is an integer division, truncated downwards

  print*,'number of blocks', nb

  i = 1    ! our current position in A
  
  do j = 1,nb-1
     call lublas2(A(i,i),bsize,n) ! LU decomposition of bsize x bsize 
                                  ! block with leading entry A(i,i)  
                                  ! Note that this block is overwritten 
                                  ! with the LU factors 
    call dtrsm('Left','Lower','No Transpose','Unit Diagonal', &
               bsize,n-j*bsize,1.0_8, &
               A(i,i),n,A(i,i+bsize),n)    ! Step 2 of the algorithm: 
                                           ! a unit lower triangular solve 
                                           ! with many right-hand sides

    call dtrsm('Right','Upper','No Transpose','Non Unit Diagonal', &
               n-j*bsize,bsize,1.0_8, &
               A(i,i),n,A(i+bsize,i),n)    ! Step 3 of the algorithm: 
                                           ! an upper  triangular solve 
                                           ! with many right-hand sides

    call dgemm('No transpose','No transpose',  &
               n-j*bsize,n-j*bsize, bsize, & 
               -1.0_8,A(i+bsize,i),n,A(i,i+bsize),n, &
                1.0_8,A(i+bsize,i+bsize),n)  ! Step 4 of the algorithm:
                                             ! form the remainder matrix 
                                             ! which is still to be factorised. 

    i = i + bsize                            ! move our position to next block
  end do 

  ! LU decomposition of final block of any remaining size
  call lublas2(A(i,i),n-i+1,n)
 
end subroutine lublas3
