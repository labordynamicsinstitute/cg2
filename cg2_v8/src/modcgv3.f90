
      SUBROUTINE MODCG( N, B, X, R, P, Q, MAXIT, RESID)
!
!     Based on CG in Dongarra "Solving Linear Systems on Vector and Shared
!      Memory Computers" p 146
      INTEGER            N,  MAXIT
      DOUBLE PRECISION   RESID
      DOUBLE PRECISION   X(N), B(N), R(N), P(N), Q(N) 


!  Purpose
!  =======
!
!  CG solves the linear system Ax = b using the
!  Conjugate Gradient iterative method with preconditioning.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER.
!          On entry, the dimension of the matrix.
!          Unchanged on exit.
!
!  B       (input) DOUBLE PRECISION array, dimension N.
!          On entry, right hand side vector B.
!          Unchanged on exit.
!
!  X       (input/output) DOUBLE PRECISION array, dimension N.
!          On input, the initial guess. This is commonly set to
!          the zero vector.
!          On exit, if INFO = 0, the iterated approximate solution.
!
!  R,W,P,Q Workspace for residual, direction vector, etc. (W not used)
!          
!
!  
!  MAXIT    (input/output) INTEGER
!          On input, the maximum iterations to be performed.
!          On output, actual number of iterations performed.
!
!  RESID   (input/output) DOUBLE PRECISION
!          On input, the allowable convergence measure for
!          norm( b - A*x ) / norm( b ).
!          On output, the final value of this measure.


!     .. Local Scalars ..
      INTEGER          ITER, ITMAX, LOGUNIT
      DOUBLE PRECISION ALPHA, BETA, RHO, RHO1, RNRM2, BNRM2, TOL, EPS
!
      EXTERNAL         DAXPY, DDOT, DNRM2
!
      LOGUNIT = 12
      EPS = 1.0E-16
      INFO = 0
      TOL = RESID
      BETA = 0.0D0

      R = B
      P = 0.0D0
!      BNRM2 = DNRM2(N,B,1)
      BNRM2 = SQRT(DOT_PRODUCT(B,B))
      IF (BNRM2 .EQ. 0) BNRM2 = 1.0
!     R <- B - AX   initial residual
      CALL MATVEC( X, R )
      R = B - R
      ITMAX = MAXIT
      MAXIT = 0
      RNRM2 =  SQRT(DOT_PRODUCT(R,R))
      RESID =  RNRM2/BNRM2
      WRITE(LOGUNIT,10) MAXIT, RNRM2, RESID
10    FORMAT(" Iteration ",I4," Norm of residual ", E15.9, " Relative error ", E15.9)
      IF ( RESID.LT.TOL ) RETURN
!      W = R
!      RHO = DOT_PRODUCT(R,W)
      RHO = DOT_PRODUCT(R,R)
!     Perform Preconditioned Conjugate Gradient iteration.
      DO MAXIT = 1, ITMAX

         P = R + BETA*P
!        Q <- A*P
         CALL MATVEC( P, Q )
         ALPHA =  RHO / DOT_PRODUCT(P,Q)
         X = X + ALPHA*P
         R = R - ALPHA*Q
         RNRM2 =  SQRT(DOT_PRODUCT(R,R))
         RESID =  RNRM2/BNRM2
         WRITE(LOGUNIT,10) MAXIT, RNRM2, RESID
         IF( RHO .LT. N*EPS) RETURN
!         IF ( RESID.LE.TOL ) RETURN
         RHO1 = RHO
!        Preconditioner Solve - Not used for Now, Just copy W <- R
!********CALL PSOLVE( W, R )
!         W = R
!         RHO = DOT_PRODUCT(R,W)
         RHO = DOT_PRODUCT(R,R)
         BETA = RHO / RHO1
      END DO
      RETURN
      END


