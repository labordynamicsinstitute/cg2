      program cg2_4v5
      use MY_DATA
      implicit none
! Revision 6/8/2006 RHC - differs from cg2_4v3 by transposing xd, xf and reading large blocks of input
      INTEGER*4 i,j,jp,jf,k,maxit,icell,info, imem,idmem,err,LRECL, nrecsperblock, iblock, itime
      INTEGER*4 persi,persim1,firmi,firmim1,logunit, int4, nvar, nmax, ierr
      INTEGER*4 flush
      REAL*8 alpha,beta,resid,yind,ybar,xsum, time(20), dclock
      REAL*4, ALLOCATABLE :: varblock(:,:)
      REAL*8, ALLOCATABLE :: covblock(:,:)
      REAL*4 yin
      EXTERNAL dclock
      CHARACTER*1 transa, transb
      CHARACTER*64 fnbase,fncgin,fncglog, fnbin4,fnbetas, fnmeans, dbgin
      itime = 1
      time(itime) = dclock()
      itime = itime + 1
      debug = .false.
      timing = .true.
      call getarg(1,fnbase)
      call getarg(2,dbgin)
      if(len(trim(dbgin)) .gt. 0) debug = .true.
      if(len(trim(fnbase)) .eq. 0) fnbase = 'cg'
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
      write(logunit,15) fncgin,fncglog,fnbin4,fnbetas,fnmeans
15    format("Input parameters file: ",A64,/, &
             "Output log file:       ",A64,/, &
             "Binary data file:      ",A64,/, &
             "Output beta file:      ",A64,/, &
             "Output means file:     ",A64,//)
      nvar = ncov + 3
      LRECL = ncov + 3                   ! recordsize in 32 bit words
      nrecsperblock = 1024*128
      iblock = 4*LRECL*nrecsperblock ! blocksize in bytes, must be less than 2GB
      ncoef = ncov + npers + nfirm - 1
      ncoef2 = ncov + nfirm -1
      ALLOCATE( dfp(ncells),dff(ncells),stat = err)
      imem = 4.0*( 2.0*dble(ncells))/(1024.0*1024.0)
      if (err .eq. 0) then
          write(logunit,20) imem
20    format("Allocated integer variables ", I10, " megabytes")
      else
         write(logunit,25) imem,err
25    format("Unable to allocate integer variables ", I10, " megabytes - error status code ",I6)
         stop
      end if
      allocate(varblock(nvar,nrecsperblock),covblock(ncov,nrecsperblock), covin(ncov),covin4(ncov), &
      covmean(ncov),d(npers),f(nfirm),df(ncells),xx(ncov,ncov), rq(ncov,ncov),xd(ncov,npers),xf(ncov,nfirm), &
           theta(ncoef),theta2(ncoef2),b(ncoef),b2(ncoef2), &
           r(ncoef2),p(ncoef2),q(ncoef2),u1(npers),stat = err)
       idmem = 8*(dble(nvar)*dble(nrecsperblock)+ ncov + 2*npers + nfirm + ncells + 2*ncov*ncov + &
              dble(npers)*dble(ncov) + dble(nfirm)*dble(ncov) + &
              5*ncoef2 +2*ncoef )/(1024.0*1024.0)
         if (err .eq. 0) then
            write(logunit,30) idmem
30          format("Allocated double variables  ", I10, " megabytes")
         else
            write(logunit,35) idmem, err
35          format("Unable to allocate double variables  ", I10, " megabytes - error status code",I6)
            stop
         end if
       idmem = idmem + imem
       write(logunit,40) idmem
40     format("Total memory allocated      ",I10, " megabytes",/)
      write(logunit,45)
45   format("Starting to read and preprocess data")
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
      time(itime) = dclock()
      write(logunit,50)  time(itime)-time(itime-1)
50    format("Finished allocating and initializing variables, times: ",2F12.2)
      itime = itime + 1
      open(unit=8,file=fnbin4,form="BINARY",recl=LRECL,status="OLD",READONLY)
      persim1 = 0
      firmim1 = 0
      ybar = 0.0D0
      covmean = 0.0D0
      k = nrecsperblock + 1    ! force reading of first block of data
      do i = 1,n
! ASSUMPTION  -  data is sorted by person, then firm
!         each record contains  yin,covin4,persi,firmi
         if (k .gt. nrecsperblock) then
            nmax = MIN(nrecsperblock,((n-i)+1))
            if (n-i-1 .le. nrecsperblock) then
               write(logunit,*) 'Reading last block - # records left = ', nmax
            end if
            read(8) varblock(:,1:nmax)
            covblock = dble(varblock(2:(ncov+1),:))
            k = 1
         end if
         covin = varblock(2:(ncov+1),k)
         yind = varblock(1,k)
         ybar = ybar + yind
         persi = transfer(varblock(ncov+2,k),int4)
         firmi = transfer(varblock(ncov+3,k),int4)
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
         xd(:,jp) = xd(:,jp) + covin
         xf(:,jf) = xf(:,jf) + covin
!         do j = 1, ncov
!            xx(:,j) = xx(:,j) + covin(j)*covin
!         end do
         if (k .eq. 1) then
            beta = 1.0D0
            call dgemm('N','T', ncov, ncov, nmax, alpha, covblock, ncov, covblock, ncov,beta, xx, ncov)
         end if
         b(1:ncov) = b(1:ncov) +yind*covin           ! accumulate X'y
         b(jp+ncov) = b(jp+ncov) + yind
         if(jf .ne. nfirm) then
            b(jf+ncov+npers) =  b(jf+ncov+npers) + yind
         end if
         persim1 = persi
         firmim1 = firmi
         k = k + 1
      end do
      close(unit=8)
      time(itime) = dclock()
      write(logunit,52) time(itime)-time(1), time(itime)-time(itime-1)
52    format("Finished reading data and filling arrays, times: ",2F12.2)
      itime = itime + 1
      ybar = ybar/dble(n)
      covmean = covmean/dble(n)
      open(unit=13,file=fnmeans,status="NEW")
      write(13,55) ybar,(covmean(i),i=1,ncov)
55    format(e15.7)
      close(unit=13)
      if (icell .ne. ncells) then
         write(logunit,60) icell,ncells
60       format(" Error - number of cells read ",I9, " not equal number specified ",I9)
      end if
     write(logunit,65)
65   format("Finished reading and preprocessing - starting preconditioning",/)
      write(logunit,70)
      if (debug) call vdmean(d,npers,xsum,          "d         ")
      if (debug) call vdmean(f,nfirm-1,xsum,        "f         ")
      if (debug) call vdmean(xd,npers*ncov,xsum,    "xd        ")
      if (debug) call vdmean(xf,nfirm*ncov,xsum,    "xf        ")
      if (debug) call vdmean(df,ncells,xsum,        "df        ")
      if (debug) call vimean(dfp,ncells,xsum,       "dfp       ")
      if (debug) call vimean(dff,ncells,xsum,       "dff       ")
      if (debug) call vdmean(b,ncoef,xsum,          "b         ")   
70    format("Means - dependent variable, covariates",/)
      write(logunit,75)  ybar,(covmean(i),i=1,ncov)
75    format(8e15.7)
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
      if (info .eq. 0) write(logunit,80)
80   format(/," Finished Cholesky part 1 DPOTRF")
      if (info .ne. 0) then
         write(logunit,85) info
85      format(" Error in DPOTRF - info = ",I5)
         stop
      end if
!    2. Save R (upper triangular part) to restore covariate effects
      rq = xx
!     3. The effect of transforming cov with inverse of u cov <- cov*inverse(u) so cov'cov = I
!         is to transform xd and xf  xd <- xd*inverse(u) xf <- xf*inverse(u)
      call dtrsm('L','U','T','N',ncov,npers,alpha,rq,ncov,xd,ncov)
      call dtrsm('L','U','T','N',ncov,nfirm,alpha,rq,ncov,xf,ncov)
!     Use diag(D'D) and diag(F'F) in preconditioning
      d = 1/sqrt(d)
      f = 1/sqrt(f)
! Commented out code from v3 when xd, xf were transposed - redo to minimize memory loads
!      do i = 1,ncov
!         xd(i,:) = xd(i,:)*d
!         xf(i,:) = xf(i,:)*f
!      end do
      do i = 1,npers
         xd(:,i) =xd(:,i)*d(i)
      end do
      do i = 1,nfirm
         xf(:,i) =xf(:,i)*f(i)
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
      call dgemm('N','T',ncov,ncov,npers,alpha,xd,ncov,xd,ncov,beta,xx,ncov)
      if (debug) call vdmean(xx,ncov*ncov,xsum,    "xx        ")
!     OLD before switching xf indicies Compute F'X-F'DD'X, store in xf
!     Compute X'F-X'DD'F, store in xf
      do i = 1,ncells
        jp = dfp(i)
        jf = dff(i)
        if (jf .le. ncoef2) then
           xf(:,jf) = xf(:,jf)-xd(:,jp)*df(i)
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
     call dgemv('N',ncov,npers,alpha,xd,ncov,b(ncov+1),1,beta,b2(1),1)
     write(logunit,*) "Finished transforming b"
      if (debug) call vdmean(b2,ncoef2,xsum,       "b2        ")      
      ierr = flush(logunit)
!    And F'Y - F'DD'Y ;D'Y is in b(ncov+1:ncov+npers), F'Y is in b2(ncov+1:ncoef2)
     do i = 1,ncells
        jp = dfp(i)+ncov
        jf = dff(i)+ncov
        if (jf .le. ncoef2) then
           b2(jf) = b2(jf) - df(i)*b(jp)
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
       time(itime) = dclock()
      write(logunit,90)  time(itime)-time(1), time(itime)-time(itime-1)     
90    format("Finished allocating and initializing variables, times: ",2F12.2)
      itime = itime + 1
       write(logunit,95)
95     format("Beginning conjugate gradient iterations")
       ierr = flush(logunit)
      call modcg(ncoef2,b2,theta2,r,p,q,maxit,resid)
      time(itime) = dclock()
      write(logunit,100)  time(itime)-time(1), time(itime)-time(itime-1)
100    format("Finished iterations , times: ",2F12.2)
      itime = itime + 1
!     Compute person effects from firm and fixed effects
!     theta = D'y - D'F psi - D'X beta
!     First D'y
      theta(ncov+1:ncov+npers) = b(ncov+1:ncov+npers)
!     Now D'y - D'Fpsi
      do i = 1,ncells
            jp = dfp(i)+ncov
            jf = dff(i)
            jf = jf+ncov
            if (jf .le. ncoef2) then
               theta(jp) = theta(jp) - theta2(jf)*df(i)
            end if
         end do   
!     Then D'y-D'Fpsi-D'Xbeta
      alpha = -1.0D0
      beta = 1.0D0
      call dgemv('T', ncov,npers, alpha, xd, ncov, theta2(1), 1, &
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
      time(itime) = dclock()
      write(logunit,105)  time(itime)-time(1), time(itime)-time(itime-1)
105   format("Finished back transformations , times: ", 2f12.2)
      itime = itime + 1
       open(unit=10,file=fnbetas,status="NEW")
      do i = 1,ncov
         j = i
         write(10,110) j,theta(i)
      end do
      do i = ncov+1, ncov+npers
         j = i-ncov
         write(10,110) j,theta(i)
      end do
      do i = ncov+npers+1,ncoef
        j = i - ncov - npers
        write(10,110) j,theta(i)
      end do
      close(unit=10)
      deallocate(dfp,dff,varblock,covin,covin4,covmean,d,f,df,xx, rq,xd,xf, &
           theta,theta2,b,b2,r,p,q,u1)
110   format(I10,1x,e15.7)
      time(itime) = dclock()
      write(logunit,120) time(itime)-time(1), time(itime)-time(itime-1)
120   format("End program , time: ",2f12.2)
      itime = itime + 1
!      stop
      end




