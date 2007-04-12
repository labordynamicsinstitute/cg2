      program cg2_4v6emix
      use MY_DATA
      
      IMPLICIT NONE
      INTEGER*4 i,j,k,jp,jf,maxit,maxvit,ivit, icell,info,imem,idmem,err,LRECL
      INTEGER*4 nrecsperblock,iblock,itime,jpers,jfirm
      INTEGER*4 persi,persim1,firmi,firmim1,logunit, int4, nvar, nvarin, nmax, ierr, flush
      INTEGER*4 invars(512)
      REAL*8 alpha,beta,resid,yind,ybar,xsum, time(20), dclock, sig2ehat, sig2ehatcg, &
           sig2phat, sig2fhat, sse, sse2, ssp, ssf
      REAL*8  degfe, degfp, degfpcg, degff, degffcg, degfecg, degfein, &
         degfpmin, degffmin, degfpmax, degffmax
      REAL*8  sig2phatlow, sig2fhatlow, sig2phathigh, sig2fhathigh, sig2pave, sig2fave, &
           lampout, lamfout, lampold, lamfold, lameps, one, zero

      CHARACTER*1 transa, transb
      REAL*4 yin
      CHARACTER*64 fnbase,fncgin,fncglog, fnbin4,fnbetas, fnmeans, fnvars, dbgin
      CHARACTER*20 varnames(512)
      CHARACTER*3  varsuffix
      LOGICAL selectvars
      REAL*4, ALLOCATABLE :: varblock(:,:)
      REAL*8, ALLOCATABLE :: covblock(:,:)
      REAL*8, ALLOCATABLE :: u1pers(:), u1ranp(:), u1firm(:), u1ranf(:), thetaranp(:), thetaranf(:)
      REAL*8, ALLOCATABLE :: b2pers(:), b2firm(:), dold(:) , fold(:)
!      lameps = .00001
! lameps is used to control convergence criterion for the lambdas.
!  since the lambdas are estimated with a Monte Carlo method, there is inherent error
!  ideally, an estimate of the error should be determined and lameps set accordingly
      lameps = .00001
      one = 1.0D0
      zero = 0.0D0
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
      fnvars = trim(fnbase) // '.vars'
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
      write(logunit,15) fncgin,fncglog,fnbin4,fnbetas,fnmeans
15    format("Input parameters file: ",A64,/, &
             "Output log file:       ",A64,/, &
             "Binary data file:      ",A64,/, &
             "Output beta file:      ",A64,/, &
             "Output means file:     ",A64,//)
      nvar = ncov + 3
      LRECL = ncov + 3                   ! recordsize in 32 bit words
      write(logunit,*) "Checking for variable selection file " // fnvars
      open(unit=7,file=fnvars,status='OLD',IOSTAT=ierr)      
      if (ierr .ne. 0) then
         selectvars = .false.
         write(logunit,*) "Using default variable names"
         write(logunit,*) "Variable name   Position in File"
         i = 1
         varnames(1) = 'Y'
         invars(1) = 1
         write(logunit,17) varnames(i), invars(i)
17       format(A20,I3)
         do i=2,ncov+1
            varnames(i) = 'X' 
            encode(3,18,varsuffix) i-1
18          format(i3)
            varnames(i) = trim(varnames(i)) // adjustl(varsuffix)
            invars(i) = i
            write(logunit,17) varnames(i), invars(i)
         end do
         i = ncov+2
         varnames(i) = 'person'
         invars(i) = i
         write(logunit,17) varnames(i), invars(i)
         i = ncov+3
         varnames(i) = 'firm'
         invars(i) = i
         write(logunit,17) varnames(i), invars(i)
         else
            write(logunit,*) "Using variables and names from "//fnvars
            write(logunit,*) "Variable name   Position in File"
         do i=1,nvar
            read(7,*) varnames(i), invars(i)
            varnames(i) = trim(varnames(i))
            write(logunit,17) varnames(i), invars(i)
         end do
! Assume firm variable is always the last variable on the record.
! If the number of variables used in this run is the same as the firm index,
! then all variables on the record are being used and we don't need to select
! variables from the input record.
         selectvars = .true.
         if (invars(nvar) .eq. nvar) selectvars = .false.
      end if
      if (selectvars) then
         nvarin = invars(nvar)
      else
         nvarin = nvar
      end if
      nrecsperblock = 1024*128
      iblock = 4*LRECL*nrecsperblock ! blocksize in bytes, must be less than 2GB
      ncoef = ncov + npers + nfirm
      ncoef2 = ncov + nfirm
      ALLOCATE(pers(n), firm(n),dfp(ncells),dff(ncells),stat = err)
      imem = 4.0*(n+n + 2.0*dble(ncells))/(1024.0*1024.0)
      if (err .eq. 0) write(logunit,20) imem
20    format("Allocated integer variables ", I10, " megabytes")
      allocate(varblock(nvarin,nrecsperblock),covblock(ncov,nrecsperblock),covin(ncov),covin4(ncov), &
           covmean(ncov),y(n),d(npers),dold(npers),din(npers),f(nfirm), fold(nfirm), fin(nfirm), &
           df(ncells),xx(ncov,ncov), xxin(ncov,ncov),rq(ncov,ncov),xd(ncov,npers),xf(ncov,nfirm), &
           theta(ncoef),theta2(ncoef2),thetaranf(ncoef2),thetaranp(ncoef2),b(ncoef),b2(ncoef2), &
           b2pers(ncoef2), b2firm(ncoef2), &
           r(ncoef2),w(ncoef2),p(ncoef2),q(ncoef2),u1(npers),u1pers(npers),u1ranp(npers), &
           u1firm(nfirm), u1ranf(nfirm), stat = err)
       idmem = 8*( dble(nvar)*dble(nrecsperblock)+ ncov + n + 2*npers + 2*nfirm + ncells + 2*ncov*ncov + &
              dble(npers)*dble(ncov) + dble(nfirm)*dble(ncov) + &
              6*ncoef +2*ncoef2 )/(1024.0*1024.0)
         if (err .eq. 0) write(logunit,30) idmem
30    format("Allocated double variables  ", I10, " megabytes")
       idmem = idmem + imem
       write(logunit,35) idmem
35     format("Total memory allocated      ",I10, " megabytes")
       thetaranp = 0.0D0
       thetaranf = 0.0D0
       maxvit = 30
       ivit = 1
      write(logunit,40)
 40   format("Starting to read and preprocess data")
      alpha = 1.0D0
      beta= 0.0D0
      icell = 0
      y = 0.0D0
      b = 0.0D0
      d = 0.0D0
      dold = 1.0D0
      din = 0.0D0
      f = 0.0D0
      fold = 1.0D0
      fin = 0.0D0
      df = 0.0D0
      xd = 0.0D0
      xf = 0.0D0
      xx = 0.0D0
      yy = 0.0D0
      time(itime) = dclock()
      write(logunit,50)  time(itime)-time(itime-1)
50    format("Finished initializing variables, times: ",2F12.2)
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
            if (selectvars) then
! Assume invars are sorted in ascending order so we don't clobber good data
               do j = 1, nvar
                  varblock(j,:) = varblock(invars(j),:)
               end do
            end if
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
      if (ivit .eq. 1) then
         open(unit=13,file=fnmeans,status="NEW")
         write(13,55) ybar,varnames(1)
         do i = 1,ncov
            write(13,55) covmean(i),varnames(i+1)
55          format(e15.7,1x,A20)
         end do
         close(unit=13)
      end if
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
      if (ivit .eq. 1) then
         write(logunit,70)
70       format("Means - dependent variable, covariates",/)      
      write(logunit,55)  ybar,varnames(1)
      do i=1,ncov
         write(logunit,55) covmean(i),varnames(i+1)
      end do
      theta = 0.0
      theta2 = 0.0
      end if
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
! Save x'x for later
      xxin = xx
      alpha = 1.0D0
      beta = 0.0D0
!     2. Compute Cholesky decomposition of xx = u'u (stored in xx)     
      call dpotrf('U',ncov,xx,ncov,info)
      if (info .eq. 0) write(logunit,14)
 14   format(" Finished Cholesky part 1 DPOTRF")
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
! NEXT ITERATION STARTS HERE - need to undo/redo preconditioning
56    d = 1/sqrt(din+lamp)
      f = 1/sqrt(fin+lamf)
      do i = 1,npers
         xd(:,i) =xd(:,i)*d(i)/dold(i)
      end do
      do i = 1,nfirm
         xf(:,i) =xf(:,i)*f(i)/fold(i)
      end do
      if (debug) call vdmean(d,npers,xsum,        "d         ")
      if (debug) call vdmean(f,nfirm,xsum,        "f         ")
      if (debug) call vdmean(xd,npers*ncov,xsum,  "xd        ")
      if (debug) call vdmean(xf,nfirm*ncov,xsum,  "xf        ")
      do i = 1,ncells
         df(i) = df(i)*d(dfp(i))*f(dff(i))/(dold(dfp(i))*fold(dff(i)))
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
      b(ncov+1:ncov+npers) =  b(ncov+1:ncov+npers)*d/dold
      b(ncov+npers+1:ncoef) = b(ncov+npers+1:ncoef)*f/fold
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
!       Transform theta to new scale 
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
       time(itime) = dclock()
      write(logunit,90)  time(itime)-time(1), time(itime)-time(itime-1)     
90    format("Finished allocating and initializing variables, times: ",2F12.2)
      itime = itime + 1
       write(logunit,19) ivit
19     format("Beginning conjugate gradient iterations",I4)
       maxit = 2000
       resid = 1.0E-10
!      call modcg(ncoef2,b2,theta2,r,w,p,q,maxit,resid)
! modcgv3.f90 does not use working vector w like modcg.f90 does
! use w for now to put in random vector to estimate degerees of freedom     
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
       degfp = npers - sig2phat/sig2p
        do jf = 1,nfirm
          j = ncov + npers + jf
          theta(j) = theta2(j-npers)*f(jf)
          sig2fhat = sig2fhat + theta(j)*theta(j)
       end do  
! Try to estimate person and firm degrees of freedom
       b2pers = 0.0D0
       if (ivit .eq. 1) then
          u1ranp = 0.0D0
          call genranvec(u1ranp,npers)
       end if
!    Now X'Y - X'D D'Y ; D'Y is in b(ncov+1:ncov+npers), X'Y is in b2(1:ncov)
     alpha = -1.0D0
     beta = 1.0D0
     call dgemv('N',ncov,npers,alpha,xd,ncov,u1ranp,1,beta,b2pers(1),1)
!    And F'Y - F'DD'Y ;D'Y is in b(ncov+1:ncov+npers), F'Y is in b2(ncov+1:ncoef2)
     do i = 1,ncells
        jpers = dfp(i)
        jfirm = dff(i)+ncov
        b2pers(jfirm) = b2pers(jfirm) - df(i)*u1ranp(jpers)
     end do
      write(logunit,*) "Beginning conjugate gradient iterations for person degrees of freedom"
       maxit = 2000
       resid = 1.0E-10
     call modcg(ncoef2,b2pers,thetaranp,r,p,q,maxit,resid)
!     Compute person effects from firm and fixed effects
!     u1pers = u1ranp - D'F psi - D'X beta
!     First D'y
      u1pers = u1ranp
!     Now D'y - D'Fpsi
      do i = 1,ncells
            jpers = dfp(i)
            jf = dff(i)
            jfirm = jf+ncov
            u1pers(jpers) = u1pers(jpers) - thetaranp(jfirm)*df(i)
         end do   
!     Then D'y-D'Fpsi-D'Xbeta
      alpha = -1.0D0
      beta = 1.0D0
      call dgemv('T', ncov,npers, alpha, xd, ncov, thetaranp(1), 1, &
              beta, u1pers, 1)
      degfpcg = dot_product(u1pers*d,d*u1ranp)* lamp
! for mixed_testdat      degfpcg = 3594.75
      ssp = sig2phat
!      sig2phat =  (ssp+sig2p*degfpcg)/npers
      sig2phat =  ssp/(npers-degfpcg)
      b2firm = 0.0D0
      if (ivit .eq. 1 ) then
         call genranvec(u1ranf,nfirm)
      end if
      b2firm((ncov+1):ncoef2) = u1ranf
            write(logunit,*) "Beginning conjugate gradient iterations for firm degrees of freedom"
       maxit = 2000
       resid = 1.0E-10
     call modcg(ncoef2,b2firm,thetaranf,r,p,q,maxit,resid)
     degffcg = dot_product(thetaranf((ncov+1):ncoef2)*f,f*u1ranf)* lamf
! for mixed_testdat     degffcg = 5128.72
      ssf  = sig2fhat
!       sig2fhat = (ssf+sig2f*degffcg)/nfirm
       sig2fhat = ssf/(nfirm-degffcg)
       degff = nfirm - sig2fhat/sig2f
!     Compute sig2ehat as Y'(Y-Xbeta-Dtheta-Fpsi)/(n-ncov)
!      and sse as (Y-Xbeta-Dtheta-Fpsi)'(Y-Xbeta-Dtheta-Fpsi) 
      sig2ehat = yy
      sse2 = yy
      alpha = -1.0D0
      beta = 1.0D0
!     now -Y'Xbeta
      do i = 1,ncov
         sig2ehat = sig2ehat - theta(i)*b(i)
         sse2 = sse2 - 2.0*theta(i)*b(i)
      end do
!     and -Y'Dtheta
      do i = 1,npers 
         sig2ehat = sig2ehat - theta(ncov+i)*b(ncov+i)/d(i)
         sse2 = sse2 - 2.0*theta(ncov+i)*b(ncov+i)/d(i)
      end do
!     and -Y'Fpsi
      do i = 1,nfirm
         sig2ehat = sig2ehat - theta(ncov+npers+i)*b(ncov+npers+i)/f(i)
         sse2 = sse2 - 2.0*theta(ncov+npers+i)*b(ncov+npers+i)/f(i)
      end do
      sig2ehat = sig2ehat/(n-ncov)
!     and beta X'X beta
      sse2 = sse2 + dot_product(theta2(1:ncov),theta2(1:ncov))
!     and theta D'D theta
      do i = 1,npers
         sse2 = sse2 + din(i)*theta(ncov+i)*theta(ncov+i)
      end do
!     and psi F'F psi
      do i = 1,nfirm
         sse2 = sse2 + fin(i)*theta(ncov+npers+i)*theta(ncov+npers+i)
      end do
!     and beta X'D theta - be careful to undo preconditioning of xd
      b2(1:ncov) = zero
      u1pers = theta(ncov+1:ncov+npers)/d
      call dgemv('N',ncov,npers,one,xd,ncov,u1pers,1,zero,b2(1),1)
      sse2 = sse2 + 2.0*dot_product(b2(1:ncov),theta2(1:ncov))
!     and beta X'F psi
      b2(1:ncov) = zero
      u1firm = theta(ncov+npers+1:ncoef)/f
! unfortunately xf contains X'F-X'DD'X not just X'F, so adjust by restoring X'F
!     Compute X'F + F'DD'X, store in xf
       do i = 1,ncells
        jp = dfp(i)
        jf = dff(i)
        xf(:,jf) = xf(:,jf)+ xd(:,jp)*df(i)
       end do 
      call dgemv('N',ncov,nfirm,one,xf,ncov,u1firm,1,zero,b2(1),1)
      sse2 = sse2 + 2.0* dot_product(b2(1:ncov),theta2(1:ncov))
!     and finally theta D'F psi
      do i = 1,ncells
        jp = dfp(i) 
        jpers = dfp(i)+ncov
        jf = dff(i)
        jfirm = dff(i)+npers+ncov
        sse2 = sse2 + 2.0*(theta(jpers)/d(jp))*df(i)*(theta(jfirm)/f(jf))
!        sse2 = sse2 + 2.0*u1pers(jp)*df(i)*u1firm(jf)
      end do          
      degfecg = ncoef - ncov - degfpcg - degffcg 
      sig2ehatcg = sse2/(n-degfecg)
       write(logunit,*) "CG estimate of error, person,firm degrees of freedom", degfecg, degfpcg,degffcg
       write(logunit,*) "CG estimate of error, person,firm variance", sig2ehatcg, sig2phat, sig2fhat    
      sse = sse2
      degfein = n - sse/sig2e
      degfe = n - sse/sig2ehat
      degfecg = ncoef - ncov - degfpcg - degffcg 
      degfpmin = 0.0D0
      do i = 1,npers
         degfpmin = degfpmin + 1.0D0/(1.0D0+din(i)/lamp)
      end do
      degffmin = 0.0D0
      do i = 1,nfirm
         degffmin = degffmin + 1.0D0/(1.0D0+fin(i)/lamf)
      end do
      lampout = sig2ehatcg/sig2phat
      lamfout = sig2ehatcg/sig2fhat
! experiments seem to show using above update formulas for lamp and lamf result
!   in fewer iterations
!     lampout = sig2ehat/sig2phat
!     lamfout = sig2ehat/sig2fhat
      write(logunit,59)  sig2ehat, sig2phat, sig2fhat, lampout, lamfout
59    format(" Updated Estimate ",/, &
             " Error variance component  ",F10.6,/, &
             " Person variance component ",F10.6,/, &
             " Firm variance component   ",F10.6,/, &
             " Person Relative Precision ",F10.6,/, &
             " Firm Relative Precision   ",F10.6,/ )
      if ((abs((lamp-lampout)/lamp) .gt. lameps) .or. (abs((lamf-lamfout)/lamf) .gt. lameps)) then
         lampold = lamp
         lamfold = lamf
         lamp = lampout
         lamf = lamfout
         dold = 1/sqrt(din+lampold)
         fold = 1/sqrt(fin+lamfold)
! restore xf
         do i =1,npers
            jp = ncov+i
            theta(jp) = theta(jp)*(din(i)+lampold)/(din(i)+lamp)
         end do
         do i =1,nfirm
            jf = ncov+npers+i
            theta(jf) = theta(jf)*(fin(i)+lamfold)/(fin(i)+lamf)
            jf = ncov+i
            thetaranf(jf) = thetaranf(jf)*f(i)*sqrt(fin(i)+lamf)
         end do
         ivit = ivit + 1
         if (ivit .le. maxvit) go to 56
      end if
      open(unit=10,file=fnbetas,status="NEW")
      do i = 1,ncov
         j = i
         write(10,60) j,theta(i), varnames(i+1)
      end do
      do i = ncov+1, ncov+npers
         j = i-ncov
         write(10,60) j,theta(i)
      end do
      do i = ncov+npers+1,ncoef
        j = i - ncov - npers
        write(10,60) j,theta(i)
      end do
 60   format(I10,1x,e15.7,1x,a20)
      time(itime) = dclock()
      write(logunit,120) time(itime)-time(1), time(itime)-time(itime-1)
120   format("End program , time: ",2f12.2)
      stop
      end
