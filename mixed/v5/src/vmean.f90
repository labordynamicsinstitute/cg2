       subroutine vdmean(x,n,xsum,xname)
         implicit none
         real*8 x(n),  xsum
         integer*4 n,  i
         character*10 xname
         xsum = sum(x)
         xsum = xsum/n
         write(12,*) "Array " // xname // " mean, n: ",xsum,n
         return
         end

       subroutine vimean(x,n,xsum,xname)
         implicit none
         real*8   xsum
         integer*4 n, x(n), i
         character*10 xname
         xsum = 0.0D0
         do i = 1, n
            xsum = xsum + dble(x(i))
         end do
         xsum = xsum/n
         write(12,*) "Array " // xname // " mean, n: ",xsum,n
         return
         end
