         subroutine MATVEC( xin,rout) 
         use my_data
         
         REAL*8 alpha, beta, xin(ncoef2),rout(ncoef2), xsum
         INTEGER*4 i,jpers,jfirm

         ! Compute the matrix vector product rout <- A*xin 
         ! with original A = (X D F)'(X D F) stored in parts in my_data
         ! In cg2, person effects D are absorbed, so the matrix A is
         ! A = (X'X-X'DD'X    X'F-X'DD'F)
         !     (F'X-F'DD'X    F'F-F'DD'F)
         ! -X'DD'X  in xx(ncov,ncov), with X'X assumed identity from preconditioning
         ! X'D in xd(npers,ncov) - not used here
         ! X'F-X'DD'F in xf(nfirm,ncov) 
         ! D'F in df(ncells), person and firm indicies in dfp(ncells) and dff(ncells)
         ! D'D is identity from preconditioning dimension(npers,npers) - not used
         ! F'F is identity from preconditioning dimension(nfirm,nfirm) 
         ! The vectors x and r have two parts:
         ! The covariate effects 1:ncov
         ! The firm effects ncov+1:ncoef2 (ncoef2=ncov+nfirm-1)
         !
         ! First the covariate effects
         rout = xin        ! xx, dd and ff are identity after transformation
         if (debug) call vdmean(xin,ncoef2,xsum,"xin       ")
         alpha = 1.0D0
         beta = 1.0D0
!        call DEGEMV(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
         call dgemv('T', nfirm-1, ncov, alpha, xf, nfirm, xin(ncov+1), 1, &
              beta, rout(1), 1)
         call dgemv('N',ncov,ncov,alpha,xx,ncov,xin(1),1,beta,rout(1),1)
         if (debug) call vdmean(rout,ncov,xsum,"rout(ncov)")
         ! Then the firm effects
         call dgemv('N', nfirm-1, ncov, alpha, xf, nfirm, xin(1), 1, &
              beta, rout(ncov+1), 1)
         if (debug) call vdmean(rout(ncov+1),nfirm-1,xsum,"rout(nfirm)")
         ! And contributions from D'F and F'D to firm effects
!   -F'DD'F for firm effects
         u1 = 0.0D0
         do i = 1,ncells
            jpers = dfp(i)
            jfirm = dff(i)+ncov
            if (jfirm .le. ncoef2) then
               u1(jpers) = u1(jpers) + xin(jfirm)*df(i)
             end if
         end do  
         if (debug) call vdmean(u1,npers,xsum,"d1        ")
         do i = 1,ncells
            jpers = dfp(i)
            jfirm = dff(i)+ncov
            if (jfirm .le. ncoef2) then
              rout(jfirm) = rout(jfirm) - u1(jpers)*df(i)
            end if
         end do 
         if (debug) call vdmean(rout,ncoef2,xsum,"rout      ")
         return
         end


