      subroutine genranvec(v,n)
! Generate a random -1,1 vector with equal # of -1's and 1's
      implicit none
      integer*4 n,i,ix
      real*8 v(n),tmp
      do i=1,n/2
         v(i) = -1.0D0
         v(n-i+1) = 1.0D0
      end do
!     take care of the odd case
      i = (n+1)/2
      v(i) = -1.0D0
      do i=1,n
         call random_number(tmp)
         ix = floor(n*tmp)+1
         tmp = v(ix)
         v(ix) = v(i)
         v(i) = tmp
      end do
      return
      end
         
         
         
