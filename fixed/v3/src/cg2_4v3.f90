      program cg2_4v3
      use MY_DATA

      INTEGER*4 i,j,jp,jf,maxit,icell,info, imem,idmem,err,LRECL
      INTEGER*4 persi,persim1,firmi,firmim1,logunit
      INTEGER*4 flush
      REAL*8 alpha,beta,resid,yind,ybar,xsum
      REAL*4 yin
      CHARACTER*1 transa, transb
      CHARACTER*64 fnbase,fncgin,fncglog, fnbin4,fnbetas, fnmeans
      debug = .false.
      call getarg(1,fnbase)
      if(len(trim(fnbase)) .eq. 0) fnbase='cg'
      fnbase = trim(fnbase)
      fncgin = trim(fnbase) // '.cgin'
      fncglog = trim(fnbase) // '.cglog'
      fnbin4 = trim(fnbase) // '.bin4'
      fnbetas = trim(fnbase) // '.betas'
      fnmeans = trim(fnbase) // '.means'
      open(unit=7,file=fncgin,status='OLD')
      read(7,*) n,ncells,npers,nfirm,ncov
      logunit = 12
      open(unit=logunit, file=fncglog,status="NEW")
      write(logunit,10) n,ncells,npers,nfirm,ncov
10    format("Number of observations   ",I10,/, &
             "Number of distinct cells ",I10,/, &
             "Number of persons        ",I10,/, &
             "Number of firms          ",I10,/, &
             "Number of covariates     ",I10,//)
      write(logunit,11) fncgin,fncglog,fnbin4,fnbetas,fnmeans
11    format("Input parameters file: ",A64,/, &
             "Output log file:       ",A64,/, &
             "Binary data file:      ",A64,/, &
             "Output beta file:      ",A64,/, &
             "Output means file:     ",A64,//)
      LRECL = ncov + 3                   ! recordsize in 32 bit words
      ncoef = ncov + npers + nfirm - 1
      ncoef2 = ncov + nfirm -1
      ALLOCATE( dfp(ncells),dff(ncells),stat = err)
      imem = 4.0*( 2.0*dble(ncells))/(1024.0*1024.0)
      if (err .eq. 0) then
          write(logunit,20) imem
20    format("Allocated integer variables ", I10, " megabytes")
      else
         write(logunit,21) imem,err
21    format("Unable to allocate integer variables ", I10, " megabytes - error status code ",I6)
         stop
      end if
      allocate(covin(ncov),covin4(ncov),covmean(ncov),d(npers),f(nfirm), &
           df(ncells),xx(ncov,ncov), rq(ncov,ncov),xd(npers,ncov),xf(nfirm,ncov), &
           theta(ncoef),theta2(ncoef2),b(ncoef),b2(ncoef2), &
           r(ncoef2),p(ncoef2),q(ncoef2),u1(npers),stat = err)
       idmem = 8*(ncov + 2*npers + nfirm + ncells + 2*ncov*ncov + &
              dble(npers)*dble(ncov) + dble(nfirm)*dble(ncov) + &
              5*ncoef2 +2*ncoef )/(1024.0*1024.0)
         if (err .eq. 0) then
            write(logunit,30) idmem
30          format("Allocated double variables  ", I10, " megabytes")
         else
            write(logunit,31) idmem, err
31          format("Unable to allocate double variables  ", I10, " megabytes - error status code",I6)
            stop
         end if
       idmem = idmem + imem
       write(logunit,35) idmem
35     format("Total memory allocated      ",I10, " megabytes",/)
      write(logunit,40)
 40   format("Starting to read and preprocess data")
      ierr = flush(logunit)
      alpha = 1.0D0
      beta= 0.0D0
      icell = 0
      b = 0.0D0
      d = 0.0D0
      f = 0.0D0
      df = 0.0D0
      xd = 0.0D0
      xf = 0.0D0
      xx = 0.0D0
      open(unit=8,file=fnbin4,form="BINARY",recl=LRECL,status="OLD")
      persim1 = 0
      firmim1 = 0
      ybar = 0.0D0
      covmean = 0.0D0
      do i = 1,n
! ASSUMPTION  -  data is sorted by person, then firm
         read(8) yin,covin4,persi,firmi
         yind = yin
         ybar = ybar + yind
         covin=covin4
         covmean = covmean + covin
         if (i .eq.1 .or. persi .ne. persim1 .or. &
         (firmi .ne. firmim1 .and. persi .eq. persim1)) &
            icell = icell + 1
         jp = persi
         d(jp) = d(jp) + 1.0D0
         dfp(icell) = jp
         jf = firmi
         f(jf) = f(jf) + 1.0D0
         dff(icell) = jf
         df(icell) = df(icell) + 1.0D0
         do j = 1, ncov
            xd(jp,j) = xd(jp,j) + covin(j)
            xf(jf,j) = xf(jf,j) + covin(j)
            do k = 1, ncov
               xx(j,k) = xx(j,k) + covin(j)*covin(k)
            end do   
         end do
         b(1:ncov) = b(1:ncov) +yind*covin           ! accumulate X'y
         b(jp+ncov) = b(jp+ncov) + yind
         if(jf .ne. nfirm) then
            b(jf+ncov+npers) =  b(jf+ncov+npers) + yind
         end if
         persim1 = persi
         firmim1 = firmi
      end do
      close(unit=8)
      ybar = ybar/dble(n)
      covmean = covmean/dble(n)
      open(unit=13,file=fnmeans,status="NEW")
      write(13,41) ybar,(covmean(i),i=1,ncov)
41    format(e15.7)
      close(unit=13)
      if (icell .ne. ncells) then
         write(logunit,42) icell,ncells
42       format(" Error - number of cells read ",I9, " not equal number specified ",I9)
      end if
     write(logunit,12)
 12   format("Finished reading and preprocessing - starting preconditioning",/)
      write(logunit,43)
      if (debug) call vdmean(d,npers,xsum,          "d         ")
      if (debug) call vdmean(f,nfirm-1,xsum,        "f         ")
      if (debug) call vdmean(xd,npers*ncov,xsum,    "xd        ")
      if (debug) call vdmean(xf,nfirm*ncov,xsum,    "xf        ")
      if (debug) call vdmean(df,ncells,xsum,        "df        ")
      if (debug) call vimean(dfp,ncells,xsum,       "dfp       ")
      if (debug) call vimean(dff,ncells,xsum,       "dff       ")
      if (debug) call vdmean(b,ncoef,xsum,          "b         ")   
43    format("Means - dependent variable, covariates",/)
      write(logunit,44)  ybar,(covmean(i),i=1,ncov)
44    format(8e15.7)
! Comment out next line if coefficients are read
       theta = 0.0
       theta2 = 0.0
!     Read saved coefficients from previous run
!
!       do i = 1, ncoef
!          read(11,*) j,theta(i)
!       end do

!
!     Preconditioning steps - use diag of D'D and F'F and a transform
!       of X'X (cov'cov) so that X'X = I
!     Transform covariates using Cholesky decomposition in several steps:
!     1. Compute qtq = cov'cov (same as X'X)
!     2. Compute Cholesky decomposition of xx = u'u (stored in xx)
!     3  Transform cov with inverse of u cov <- cov*inverse(u) so cov'cov = I
!     4. Save R (upper triangular part) to restore covariate effects
!    
!     This version absorbs person effects into firm effects and
!     the preconditioning includes the absorption steps.

 

!

      alpha = 1.0D0
      beta = 0.0D0
!     2. Compute Cholesky decomposition of xx = u'u (stored in xx)     
      call dpotrf('U',ncov,xx,ncov,info)
      if (info .eq. 0) write(logunit,15)
 15   format(/," Finished Cholesky part 1 DPOTRF")
      if (info .ne. 0) then
         write(logunit,16) info
 16      format(" Error in DPOTRF - info = ",I5)
         stop
      end if
!    2. Save R (upper triangular part) to restore covariate effects
      rq = xx
!      write(12,14) xx
!     3. The effect of transforming cov with inverse of u cov <- cov*inverse(u) so cov'cov = I
!         is to transform xd and xf  xd <- xd*inverse(u) xf <- xf*inverse(u)
!      The following line was in testcg (cg) - not used in iterations but used in b2
      call dtrsm('R','U','N','N',npers,ncov,alpha,rq,ncov,xd,npers)
      call dtrsm('R','U','N','N',nfirm,ncov,alpha,rq,ncov,xf,nfirm)
!     Use diag(D'D) and diag(F'F) in preconditioning
      d = 1/sqrt(d)
      f = 1/sqrt(f)
      do i = 1,ncov
         xd(:,i) = xd(:,i)*d
         xf(:,i) = xf(:,i)*f
      end do
      if (debug) call vdmean(d,npers,xsum,        "d         ")
      if (debug) call vdmean(f,nfirm,xsum,        "f         ")
      if (debug) call vdmean(xd,npers*ncov,xsum,  "xd        ")
      if (debug) call vdmean(xf,nfirm*ncov,xsum,  "xf        ")
      do i = 1,ncells
         df(i) = df(i)*d(dfp(i))*f(dff(i))
      end do
      if (debug) call vdmean(df,ncells,xsum,      "df        ")
!     Transformed X'X is now identity
      xx = 0.0D0
!     Compute X'X-X'DD'X - X'X is Identity, store just -X'DD'X in xx
!SUBROUTINE DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,BETA, C, LDC)
      alpha = -1.0D0
      beta = 0.0D0
      call dgemm('T','N',ncov,ncov,npers,alpha,xd,npers,xd,npers,beta,xx,ncov)
      if (debug) call vdmean(xx,ncov*ncov,xsum,    "xx        ")
!     Compute F'X-F'DD'X, store in xf
      do i = 1,ncells
        jpers = dfp(i)
        jfirm = dff(i)
        if (jfirm .le. ncoef2) then
           do j = 1,ncov
              xf(jfirm,j) = xf(jfirm,j)-xd(jpers,j)*df(i)
           end do
        end if
      end do  
      if (debug) call vdmean(xf,ncov*nfirm,xsum,"xf        ")
!     b <- (X D F)'y - already done just need to condition consistently
      alpha = 1.0D0
      call dtrsm('R','U','N','N',1,ncov,alpha,rq,ncov,b,1)
      b2(1:ncov) = b(1:ncov)
      b(ncov+1:ncov+npers) =  b(ncov+1:ncov+npers)*d
      b(ncov+npers+1:ncoef) = b(ncov+npers+1:ncoef)*f(1:nfirm-1)
      b2(ncov+1:ncoef2) = b(ncov+npers+1:ncoef)
      if (debug) call vdmean(b,ncoef,xsum,         "b         ")   
      if (debug) call vdmean(b2,ncoef2,xsum,       "b2         ")   
!    Now X'Y - X'DD'Y ; D'Y is in b(ncov+1:ncov+npers), X'Y is in b2(1:ncov)
     alpha = -1.0D0
     beta = 1.0D0
     call dgemv('T',npers,ncov,alpha,xd,npers,b(ncov+1),1,beta,b2(1),1)
     write(logunit,*) "Finished transforming b"
      if (debug) call vdmean(b2,ncoef2,xsum,       "b2        ")      
      ierr = flush(logunit)
!    And F'Y - F'DD'Y ;D'Y is in b(ncov+1:ncov+npers), F'Y is in b2(ncov+1:ncoef2)
     do i = 1,ncells
        jpers = dfp(i)+ncov
        jfirm = dff(i)+ncov
        if (jfirm .le. ncoef2) then
           b2(jfirm) = b2(jfirm) - df(i)*b(jpers)
        end if
      end do
      if (debug) call vdmean(b2,ncoef2,xsum,       "b2         ")   
      resid = 1.0E-10
      maxit = 2000
!       Transform theta to new scale if old theta read
!     CALL DTRMV(UPLO, TRANS, DIAG, N, A, LDA, X, INCX)
      call dtrmv('U','N','N', ncov, rq, ncov, theta, 1)
      theta2(1:ncov) = theta(1:ncov)
       do jp = 1,npers
          j = ncov + jp
          theta(j) = theta(j)/d(jp)
       end do  
       do jf = 1,nfirm-1
          j = ncov + npers + jf
          theta(j) = theta(j)/f(jf)
          theta2(j-npers) = theta(j)
       end do
       write(logunit,17)
17     format("Beginning conjugate gradient iterations")
       ierr = flush(logunit)
      call modcg(ncoef2,b2,theta2,r,p,q,maxit,resid)
!     Compute person effects from firm and fixed effects
!     theta = D'y - D'F psi - D'X beta
!     First D'y
      theta(ncov+1:ncov+npers) = b(ncov+1:ncov+npers)
!     Now D'y - D'Fpsi
      do i = 1,ncells
            jpers = dfp(i)+ncov
            jf = dff(i)
            jfirm = jf+ncov
            if (jfirm .le. ncoef2) then
               theta(jpers) = theta(jpers) - theta2(jfirm)*df(i)
            end if
         end do   
!     Then D'y-D'Fpsi-D'Xbeta
      alpha = -1.0D0
      beta = 1.0D0
      call dgemv('N', npers, ncov, alpha, xd, npers, theta2(1), 1, &
              beta, theta(ncov+1), 1)
 !     Transform theta back to original scale
       theta(1:ncov) = theta2(1:ncov)
       call dtrsv('U','N','N', ncov, rq, ncov,theta(1),1)
        do jp = 1,npers
          j = ncov + jp
          theta(j) = theta(j)*d(jp)
       end do  
       do jf = 1,nfirm-1
          j = ncov + npers + jf
          theta(j) = theta2(j-npers)*f(jf)
       end do  
       open(unit=10,file=fnbetas,status="NEW")
      do i = 1,ncov
         j = i
         write(10,60) j,theta(i)
      end do
      do i = ncov+1, ncov+npers
         j = i-ncov
         write(10,60) j,theta(i)
      end do
      do i = ncov+npers+1,ncoef
        j = i - ncov - npers
        write(10,60) j,theta(i)
      end do
      close(unit=10)
 60   format(I10,1x,e15.7)
      stop
      end





