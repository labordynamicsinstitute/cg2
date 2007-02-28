      program cg2_4v5mix
      use MY_DATA
      
      IMPLICIT NONE
      INTEGER*4 i,j,k,jp,jf,maxit,icell,info,imem,idmem,err,LRECL,nrecsperblock,iblock,itime,jpers,jfirm
      INTEGER*4 persi,persim1,firmi,firmim1,logunit, int4, nvar, nmax, ierr, flush
      REAL*8 alpha,beta,resid,yind,ybar,xsum, time(20), dclock, sig2ehat, sig2phat, sig2fhat, sse, ssp, ssf
      REAL*8  degfe, degfp, degff, degfein, degfpmin, degffmin, degfpmax, degffmax
      REAL*8  sig2phatlow, sig2fhatlow, sig2phathigh, sig2fhathigh, sig2pave, sig2fave, lampout, lamfout
      CHARACTER*1 transa, transb
      REAL*4 yin
      CHARACTER*64 fnbase,fncgin,fncglog, fnbin4,fnbetas, fnmeans, dbgin
      REAL*4, ALLOCATABLE :: varblock(:,:)
      REAL*8, ALLOCATABLE :: covblock(:,:)
      itime = 1
      time(itime) = dclock()
      itime = itime + 1
      debug = .false.
      timing = .true.
      call getarg(1,fnbase)
      call getarg(2,dbgin)
      if(len(trim(dbgin)) .gt. 0) debug = .true.
      if(len(trim(fnbase)) .eq. 0) fnbase='cg'
      fnbase = trim(fnbase)
      fncgin = trim(fnbase) // '.cgmixin'
      fncglog = trim(fnbase) // '.cgmixlog'
      fnbin4 = trim(fnbase) // '.bin4'
      fnbetas = trim(fnbase) // '.mixbetas'
      fnmeans = trim(fnbase) // '.mixmeans'
      open(unit=7,file=fncgin,status='OLD')      
      read(7,*) n,ncells,npers,nfirm,ncov, sig2e, sig2p, sig2f
      lamp = sig2e/sig2p
      lamf = sig2e/sig2f
      logunit = 12
      open(unit=logunit, file=fncglog,status="NEW")
      write(logunit,10) n,ncells,npers,nfirm,ncov, sig2e, sig2p, sig2f, lamp, lamf
10    format("Number of observations   ",I10,/, &
             "Number of distinct cells ",I10,/, &
             "Number of persons        ",I10,/, &
             "Number of firms          ",I10,/, &
             "Number of covariates     ",I10,/, &
             " Error variance component  ",F10.6,/, &
             " Person variance component ",F10.6,/, &
             " Firm variance component   ",F10.6,/, &
             " Person Relative Precision ",F10.6,/, &
             " Firm Relative Precision   ",F10.6,/ )
      nvar = ncov + 3
      LRECL = ncov + 3                   ! recordsize in 32 bit words
      nrecsperblock = 1024*128
      iblock = 4*LRECL*nrecsperblock ! blocksize in bytes, must be less than 2GB
      ncoef = ncov + npers + nfirm
      ncoef2 = ncov + nfirm
      ALLOCATE(pers(n), firm(n),dfp(ncells),dff(ncells),stat = err)
      imem = 4.0*(n+n + 2.0*dble(ncells))/(1024.0*1024.0)
      if (err .eq. 0) write(logunit,20) imem
20    format("Allocated integer variables ", I10, " megabytes")
      allocate(varblock(nvar,nrecsperblock),covblock(ncov,nrecsperblock),covin(ncov),covin4(ncov), &
           covmean(ncov),y(n),d(npers),din(npers),f(nfirm), fin(nfirm), &
           df(ncells),xx(ncov,ncov), rq(ncov,ncov),xd(ncov,npers),xf(ncov,npers), &
           theta(ncoef),theta2(ncoef2),b(ncoef),b2(ncoef2), &
           r(ncoef2),w(ncoef2),p(ncoef2),q(ncoef2),u1(npers),stat = err)
       idmem = 8*( dble(nvar)*dble(nrecsperblock)+ ncov + n + 2*npers + 2*nfirm + ncells + 2*ncov*ncov + &
              dble(npers)*dble(ncov) + dble(nfirm)*dble(ncov) + &
              6*ncoef +2*ncoef2 )/(1024.0*1024.0)
         if (err .eq. 0) write(logunit,30) idmem
30    format("Allocated double variables  ", I10, " megabytes")
       idmem = idmem + imem
       write(logunit,35) idmem
35     format("Total memory allocated      ",I10, " megabytes")
      write(logunit,40)
 40   format("Starting to read and preprocess data")
      alpha = 1.0D0
      beta= 0.0D0
      icell = 0
      y = 0.0D0
      b = 0.0D0
      d = 0.0D0
      f = 0.0D0
      df = 0.0D0
      xd = 0.0D0
      xf = 0.0D0
      xx = 0.0D0
      yy = 0.0D0
      time(itime) = dclock()
      write(logunit,50)  time(itime)-time(itime-1)
50    format("Finished allocating and initializing variables, times: ",2F12.2)
!  Remove convert='BIG_ENDIAN' when running against files created on the current computer
!   or any other machine with the same ENDIANNESS. Also remove it in  open farther below fro reread
!      open(unit=8,file=fnbin4,form="BINARY",recl=LRECL,convert='BIG_ENDIAN',status="OLD")
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
         pers(i) = transfer(varblock(ncov+2,k),int4)
         firm(i) = transfer(varblock(ncov+3,k),int4)
         covmean = covmean + covin
         y(i) = yind
         yy = yy + y(i)*y(i)
         if (i .eq.1 .or. pers(i) .ne. pers(i-1) .or. &
         (firm(i) .ne. firm(i-1) .and. pers(i) .eq. pers(i-1))) &
            icell = icell + 1
         jp = pers(i)
         din(jp) = din(jp) + 1.0D0
         dfp(icell) = jp
         jf = firm(i)
         fin(jf) = fin(jf) + 1.0D0
         dff(icell) = jf
         df(icell) = df(icell) + 1.0D0
 
         xd(:,jp) = xd(:,jp) + covin
         xf(:,jf) = xf(:,jf) + covin
         if (k .eq. 1) then
            beta = 1.0D0
!           call dgemm('N','T', ncov, ncov, nmax, alpha, covblock, ncov, covblock, ncov,beta, xx, ncov)
            call dsyrk('U','N', ncov, nmax, alpha, covblock, ncov, beta, xx, ncov)
         end if       
         b(1:ncov) = b(1:ncov) +y(i)*covin           ! accumulate X'y
         b(jp+ncov) = b(jp+ncov) + y(i)
         b(jf+ncov+npers) =  b(jf+ncov+npers) + y(i)
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
      f = fin
      d = din
      if (icell .ne. ncells) then
         write(logunit,42) icell,ncells
42       format(" Error - number of cells read ",I9, " not equal number specified ",I9)
      end if

      if (debug) call vdmean(d,npers,xsum,          "d         ")
      if (debug) call vdmean(f,nfirm-1,xsum,        "f         ")
      if (debug) call vdmean(xd,npers*ncov,xsum,    "xd        ")
      if (debug) call vdmean(xf,nfirm*ncov,xsum,    "xf        ")
      if (debug) call vdmean(df,ncells,xsum,        "df        ")
      if (debug) call vimean(dfp,ncells,xsum,       "dfp       ")
      if (debug) call vimean(dff,ncells,xsum,       "dff       ")
      if (debug) call vdmean(b,ncoef,xsum,          "b         ")   

      write(logunit,70)
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

      write(logunit,12)
 12   format(" Finished reading and preprocessing - starting preconditioning")
      close(unit=8)
!

      alpha = 1.0D0
      beta = 0.0D0
!     2. Compute Cholesky decomposition of xx = u'u (stored in xx)     
      call dpotrf('U',ncov,xx,ncov,info)
      if (info .eq. 0) write(logunit,15)
 15   format(" Finished Cholesky part 1 DPOTRF")
      if (info .ne. 0) then
         write(logunit,16) info
 16      format(" Error in DPOTRF - info = ",I5)
         stop
      end if
!    2. Save R (upper triangular part) to restore covariate effects
      rq = xx
!     3. The effect of transforming cov with inverse of u cov <- cov*inverse(u) so cov'cov = I
!         is to transform xd and xf  xd <- xd*inverse(u) xf <- xf*inverse(u)
      call dtrsm('L','U','T','N',ncov,npers,alpha,rq,ncov,xd,ncov)
      call dtrsm('L','U','T','N',ncov,nfirm,alpha,rq,ncov,xf,ncov)
!     Use diag(D'D) and diag(F'F) in preconditioning
      d = 1/sqrt(d+lamp)
      f = 1/sqrt(f+lamf)
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
!     Compute F'X-F'DD'X, store in xf
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
      b2(1:ncov) = b(1:ncov)
      call dtrsm('R','U','N','N',1,ncov,alpha,rq,ncov,b2,1)
      b(ncov+1:ncov+npers) =  b(ncov+1:ncov+npers)*d
      b(ncov+npers+1:ncoef) = b(ncov+npers+1:ncoef)*f(1:nfirm)
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
        jpers = dfp(i)+ncov
        jfirm = dff(i)+ncov
        b2(jfirm) = b2(jfirm) - df(i)*b(jpers)
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
       do jf = 1,nfirm
          j = ncov + npers + jf
          theta(j) = theta(j)/f(jf)
          theta2(j-npers) = theta(j)
       end do

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
       do jf = 1,nfirm
          j = ncov + npers + jf
          theta(j) = theta(j)/f(jf)
          theta2(j-npers) = theta(j)
       end do
       time(itime) = dclock()
      write(logunit,90)  time(itime)-time(1), time(itime)-time(itime-1)     
90    format("Finished allocating and initializing variables, times: ",2F12.2)
      itime = itime + 1
       write(logunit,17)
17     format("Beginning conjugate gradient iterations")
!      call modcg(ncoef2,b2,theta2,r,w,p,q,maxit,resid)
! modcgv3.f90 does not use working vector w like modcg.f90 does
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
            theta(jpers) = theta(jpers) - theta2(jfirm)*df(i)
         end do   
!     Then D'y-D'Fpsi-D'Xbeta
      alpha = -1.0D0
      beta = 1.0D0
      call dgemv('T', ncov,npers, alpha, xd, ncov, theta2(1), 1, &
              beta, theta(ncov+1), 1)
 !     Transform theta back to original scale
       theta(1:ncov) = theta2(1:ncov)
       call dtrsv('U','N','N', ncov, rq, ncov,theta(1),1)
       sig2ehat = 0.0D0
       sig2phat = 0.0D0
       sig2fhat = 0.0D0
        do jp = 1,npers
          j = ncov + jp
          theta(j) = theta(j)*d(jp)
          sig2phat = sig2phat + theta(j)*theta(j)
       end do
       ssp = sig2phat
       degfp = npers - sig2phat/sig2p
       sig2phat = sig2phat/npers
       do jf = 1,nfirm
          j = ncov + npers + jf
          theta(j) = theta2(j-npers)*f(jf)
          sig2fhat = sig2fhat + theta(j)*theta(j)
       end do  
       ssf  = sig2fhat
       degff = nfirm - sig2fhat/sig2f
       sig2fhat = sig2fhat/nfirm
!     Compute sig2ehat as Y'(Y-Xbeta-Dtheta-Fpsi)/(n-ncov)
      sig2ehat = yy
      sse = 0.0D0
      alpha = -1.0D0
      beta = 1.0D0
!     now -Y'Xbeta
      do i = 1,ncov
         sig2ehat = sig2ehat - theta(i)*b(i)
      end do
!     and -Y'Dtheta
      do i = 1,npers 
         sig2ehat = sig2ehat - theta(ncov+i)*b(ncov+i)/d(i)
      end do
!     and -Y'Fpsi
      do i = 1,nfirm
         sig2ehat = sig2ehat - theta(ncov+npers+i)*b(ncov+npers+i)/f(i)
      end do
      sig2ehat = sig2ehat/(n-ncov)
      write(logunit,19)
19    format("Rereading data to compute residuals and diagnostics")
!      open(unit=8,file=fnbin4,form="BINARY",recl=LRECL,convert='BIG_ENDIAN',status="OLD")
      open(unit=8,file=fnbin4,form="BINARY",recl=LRECL,status="OLD")
      do i=1,n
         read(8) yin,covin4,pers(i),firm(i)
         y(i) = yin
         covin = covin4
         jp = pers(i)
         jf = firm(i)
         y(i) = y(i) - theta(ncov+jp) - theta(ncov+npers+jf)
         do j = 1,ncov
            y(i) = y(i) - covin(j)*theta(j)
         end do
         sse = sse + y(i)*y(i)
      end do
      degfein = n - sse/sig2e
      degfe = n - sse/sig2ehat
      degfpmin = 0.0D0
      do i = 1,npers
         degfpmin = degfpmin + 1.0D0/(1.0D0+din(i)/lamp)
      end do
      degffmin = 0.0D0
      do i = 1,nfirm
         degffmin = degffmin + 1.0D0/(1.0D0+fin(i)/lamf)
      end do
      degfpmax = (ncov + npers + nfirm) - degfe - degffmin
      degffmax = (ncov + npers + nfirm) - degfe - degfpmin
      sig2phatlow = ssp/(npers-degfpmin)
      sig2fhatlow = ssf/(nfirm-degffmin)
      sig2phathigh = ssp/(npers-degfpmax)
      sig2fhathigh = ssf/(nfirm-degffmax)
      sig2pave = .5*(sig2phatlow+sig2phathigh)
      sig2fave = .5*(sig2fhatlow+sig2fhathigh)
      lampout = sig2ehat/sig2pave
      lamfout = sig2ehat/sig2fave
      write(logunit,59) sig2e, sig2ehat, degfein,degfe, &
           sig2p, sig2phatlow, sig2pave, sig2phathigh, lamp, lampout, &
           degfp, degfpmin, .5*(degfpmin+degfpmax),degfpmax, &
           sig2f, sig2fhatlow, sig2fave, sig2fhathigh, lamf, lamfout, &
           degff, degffmin, .5*(degffmin+degffmax),degffmax, &
           sse, ssp, ssf
59    format(" Error variance component, - input, estimate                      ",2F10.6,/, &
             " Error degrees of freedom                                         ",2F14.2,/, &
             " Person variance component - input, low , average, high estimate  ",4F10.6,/, &
             " Person relative precision - input, average out                   ",2F10.6,/, &
             " Person degees of freedom -  input, low , average, high estimate  ",4F14.2,/, &
             " Firm variance component -   input, low , average, high estimate  ",4F10.6,/, &
             " Firm relative precision - input, average out                     ",2F10.6,/, &
             " Firm degees of freedom -    input, low , average, high estimate  ",4F14.2,/, &
             " Error sum of squares      ",F10.1,/, &
             " Person sum of squares     ",F10.2,/, &
             " Firm sum of squares       ",F10.2,/)
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
 60   format(I10,1x,e15.7)
      stop
      end
