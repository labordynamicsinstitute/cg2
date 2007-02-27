      program cg2_4mix
      use MY_DATA
      
      IMPLICIT NONE
      INTEGER*4 i,j,k,jp,jf,maxit,icell,info, imem,idmem,err,LRECL,jpers, jfirm, logunit
      REAL*8 alpha,beta,resid, sig2ehat, sig2phat, sig2fhat, sse, ssp, ssf
      REAL*8  degfe, degfp, degff, degfein, degfpmin, degffmin, degfpmax, degffmax
      REAL*8  sig2phatlow, sig2fhatlow, sig2phathigh, sig2fhathigh, sig2pave, sig2fave, lampout, lamfout
      CHARACTER*1 transa, transb
      REAL*4 yin
      LOGICAL*1 debug
      CHARACTER*64 fnbase,fncgin,fncglog, fnbin4,fnbetas, fnmeans
      debug = .false.
      call getarg(1,fnbase)
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
      LRECL = 2*ncov + 4                    ! recordsize in 32 bit words
      ncoef = ncov + npers + nfirm
      ncoef2 = ncov + nfirm
      ALLOCATE(pers(n), firm(n),dfp(ncells),dff(ncells),stat = err)
      imem = 4*(n+n + 2*ncells)/(1024*1024)
      if (err .eq. 0) write(logunit,20) imem
20    format("Allocated integer variables ", I10, " megabytes")
      allocate(covin(ncov),covin4(ncov),y(n),d(npers),din(npers),f(nfirm), fin(nfirm), &
           df(ncells),xx(ncov,ncov), rq(ncov,ncov),xd(npers,ncov),xf(nfirm,ncov), &
           theta(ncoef),theta2(ncoef2),b(ncoef),b2(ncoef2), &
           r(ncoef2),w(ncoef2),p(ncoef2),q(ncoef2),u1(npers),stat = err)
       idmem = 8*(ncov + n + 2*npers + nfirm + ncells + 2*ncov*ncov + &
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
!  Remove convert='BIG_ENDIAN' when running against files created on the current computer
!   or any other machine with the same ENDIANNESS. Also remove it in  open farther below fro reread
!      open(unit=8,file=fnbin4,form="BINARY",recl=LRECL,convert='BIG_ENDIAN',status="OLD")
      open(unit=8,file=fnbin4,form="BINARY",recl=LRECL,status="OLD")
      do i = 1,n
! ASSUMPTION  -  data is sorted by person, then firm
         read(8) yin,covin4,pers(i),firm(i)
         y(i) = yin
         covin = covin4
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
         do j = 1, ncov
            xd(jp,j) = xd(jp,j) + covin(j)
            xf(jf,j) = xf(jf,j) + covin(j)
            do k = 1, ncov
               xx(j,k) = xx(j,k) + covin(j)*covin(k)
            end do   
         end do
         b(1:ncov) = b(1:ncov) +y(i)*covin           ! accumulate X'y
         b(jp+ncov) = b(jp+ncov) + y(i)
         b(jf+ncov+npers) =  b(jf+ncov+npers) + y(i)
      end do
      f = fin
      d = din
      if (icell .ne. ncells) then
         write(logunit,42) icell,ncells
42       format(" Error - number of cells read ",I9, " not equal number specified ",I9)
      end if
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
!      write(12,14) xx
!     3. The effect of transforming cov with inverse of u cov <- cov*inverse(u) so cov'cov = I
!         is to transform xd and xf  xd <- xd*inverse(u) xf <- xf*inverse(u)
!      The following line was in testcg (cg) - not used in iterations but used in b2
      call dtrsm('R','U','N','N',npers,ncov,alpha,rq,ncov,xd,npers)
      call dtrsm('R','U','N','N',nfirm,ncov,alpha,rq,ncov,xf,nfirm)
!     Use diag(D'D) and diag(F'F) in preconditioning
      d = 1/sqrt(d+lamp)
      f = 1/sqrt(f+lamf)
      do i = 1,ncov
         xd(:,i) = xd(:,i)*d
         xf(:,i) = xf(:,i)*f
      end do
      do i = 1,ncells
         df(i) = df(i)*d(dfp(i))*f(dff(i))
      end do
!     Transformed X'X is now identity
      xx = 0.0D0
!     Compute X'X-X'DD'X - X'X is Identity, store just -X'DD'X in xx
!SUBROUTINE DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,BETA, C, LDC)
      alpha = -1.0D0
      beta = 0.0D0
      call dgemm('T','N',ncov,ncov,npers,alpha,xd,npers,xd,npers,beta,xx,ncov)
!     Compute F'X-F'DD'X, store in xf
      do i = 1,ncells
        jpers = dfp(i)
        jfirm = dff(i)
        do j = 1,ncov
          xf(jfirm,j) = xf(jfirm,j)-xd(jpers,j)*df(i)
        end do
      end do  
!     b <- (X D F)'y - already done just need to condition consistently
      alpha = 1.0D0
      b2(1:ncov) = b(1:ncov)
      call dtrsm('R','U','N','N',1,ncov,alpha,rq,ncov,b2,1)
      b(ncov+1:ncov+npers) =  b(ncov+1:ncov+npers)*d
      b(ncov+npers+1:ncoef) = b(ncov+npers+1:ncoef)*f(1:nfirm)
      b2(ncov+1:ncoef2) = b(ncov+npers+1:ncoef)
!    Now X'Y - X'DD'Y ; D'Y is in b(ncov+1:ncov+npers), X'Y is in b2(1:ncov)
     alpha = -1.0D0
     beta = 1.0D0
     call dgemv('T',npers,ncov,alpha,xd,npers,b(ncov+1),1,beta,b2(1),1)
!    And F'Y - F'DD'Y ;D'Y is in b(ncov+1:ncov+npers), F'Y is in b2(ncov+1:ncoef2)
     do i = 1,ncells
        jpers = dfp(i)+ncov
        jfirm = dff(i)+ncov
        b2(jfirm) = b2(jfirm) - df(i)*b(jpers)
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
      call dgemv('N', npers, ncov, alpha, xd, npers, theta2(1), 1, &
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
      write(logunit,50) sig2e, sig2ehat, degfein,degfe, &
           sig2p, sig2phatlow, sig2pave, sig2phathigh, lamp, lampout, &
           degfp, degfpmin, .5*(degfpmin+degfpmax),degfpmax, &
           sig2f, sig2fhatlow, sig2fave, sig2fhathigh, lamf, lamfout, &
           degff, degffmin, .5*(degffmin+degffmax),degffmax, &
           sse, ssp, ssf
50    format(" Error variance component, - input, estimate                      ",2F10.6,/, &
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
