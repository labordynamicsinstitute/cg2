* $Id$
* $URL$
********************************************************
*
* 	Conjugate gradient
*
* Based on work by Rob Creecy
* Stata version by Amine Ouazad
********************************************************

********************************************************
** Loading data for CG algorithm
********************************************************



do params.do

use $CGIN
sort $personid $firmid

set more off 

mata: mata clear

defvar

mata

function maincg(string rowvector covariates) {

	real scalar i, j, jp, jf;
	real vector dfp, dff, df, d;
	real matrix xd, xf, xx;
	real matrix u;	
	real vector b;

	real vector theta;

	printf("Loading data in Mata ...\n");

	printf("Covariates\n");
	st_view(cov=., 	., covariates)

	printf("Dependent variable\n");
	st_view(y=., 	., "$DEPENDENT")

	printf("Person IDs\n");
	st_view(person=., 	., "$personid")

	printf("Firm IDs\n");
	st_view(firm=., ., "$firmid")

	printf("Dependent variable : $DEPENDENT\n");

	printf("Problem Size ...\n");

	n 		= 	length(y);
	ncov		=	cols(cov);
	npersons 	=	$NPERSONS;
	nfirms	=	$NFIRMS;
	ncells	=	$NCELLS;
	ncoef		= 	npersons + nfirms + ncov - 1;

	printf("%f observations, %f covariates, %f persons, %f firms, %f cells\n",n,ncov,npersons,nfirms,ncells);

	if ((length(y) != rows(cov)) | (length(y)!=length(person)) | (length(y)!=length(firm)) ) {
		printf("Error : invalid matrix length\n");
		exit();
	}

	printf("Initializing normal form matrices\n");

	icell = 1;
	xd	= J(npersons,ncov,0);
	xf	= J(nfirms,ncov,0);
	d	= J(npersons,1, 0);
	dfp	= J(ncells, 1, 0);
	dff	= J(ncells, 1, 0);
	df	= J(ncells, 1, 0);
	f	= J(nfirms, 1, 0);

	printf("Starting loop\n");

	for ( i=1 ; i<=n ; i++ ) {

		if (i>1) {
			if (person[i] != person[i-1] || (firm[i] != firm[i-1] && person[i] == person[i-1])) {
				icell = icell +1;
			}
		}

		/* Construction de d */
		jp = person[i]
		d[jp] = d[jp] + 1

		/* Construction de dfp */
		dfp[icell] = jp

		/* Construction de f */
		jf = firm[i]
		f[jf] = f[jf] + 1

		/* Construction de dff */
		dff[icell] = jf

		/* Construction de df */
		df[icell] = df[icell] + 1

		for ( j = 1; j <= ncov ; j++ ) {
		/* Construction de XD et XF */
	
			xd[jp,j] = xd[jp, j] + cov[i, j]
			xf[jf,j] = xf[jf, j] + cov[i, j]
		}
	}

	if (icell != ncells) {
		printf("Error number of cells read %f not equal to %f\n", icell, ncells);
		exit();
	}

	theta = J(ncoef, 1, 0);

	/* 	Preconditioning steps - use diag of D'D and F'F and a transform
		of X'X (cov'cov) so that X'X = I
		Transform covariates using Cholesky Decomposition in several steps
		Compute Qtq = cov'cov
		Compute Cholesky Decomposition of xx=u'u
		Transform COV with inverse of u cov <- cov * inverse(u) so cov'cov = I
		Save R (upper triangular part) to restore covariate effects 
	*/
	
	printf("Finished preprocessing - starting preconditioning\n");

	xx = cross(cov,cov);
	printf("Finished COV'COV\n");

	xx;
	
	xx = cholesky(xx);
	printf("Finished Cholesky decomposition of COV'COV\n");
	
	u = xx';

	u;

	luinv(u);

	cov 				= cov * luinv(u) ;
	xd 				= xd * luinv(u) ;
	printf("Finished XD<-XD*u^-1\n");

	xf 				= xf * luinv(u) ;
	printf("Finished XF<-XF*u^-1\n");

	for (i = 1; i<=npersons ; i ++) {
		d[i] = 1/sqrt(d[i]);
	}

	for (i = 1; i<=nfirms ; i ++) {
		f[i] = 1/sqrt(f[i]);
	}

	for (j=1;j<=ncov;j++) {
		for (i = 1; i<=npersons ; i ++ ) {
			xd[i, j] = xd[i, j]*d[i];
		}
		for (i = 1; i<=nfirms ; i ++ ) {
			xf[i, j] = xf[i, j]*f[i];
		}
	}

	for (i=1; i<=ncells; i++) {	
		df[i] = df[i]*d[dfp[i]]*f[dff[i]]
	}

	// Does (X D F)'y and puts it in b 
	xtprod (person, firm, cov[1 .. n,1..ncov], f, d, y, b, ncov, npersons, nfirms, ncells, n);
	b = b'

	resid = 1.0e-7
	maxit = 1000

	printf("Scaling theta for preconditioning\n");

	u
	theta[1 .. ncov]
	
	theta[1..ncov] = u * theta[1..ncov]

	printf("Person effects\n");

	for (jp = 1; jp<=npersons; jp++) {
		j = ncov + jp
		theta[j]=theta[j]/d[jp]	
	}

	printf("Firm effects\n");

	for (jf = 1; jf<=nfirms-1; jf++) {
		j = ncov + npersons +jf
		theta[j] = theta[j]/f[jf]
	}
	
	printf("Beginning Conjugate Gradient Iterations\n");

	modcg(b,theta,maxit,resid, xd, xf, df, dfp, dff, ncov, npersons, nfirms, ncells, ncoef, n);

	printf("Transforming theta back to original scale\n");

	printf("Size of u (%f x %f) and of theta (%f x %f)\n", rows(u), cols(u), rows(theta), cols(theta));

	theta[1 .. ncov] = luinv(u) * theta[1 .. ncov]

	for (jp = 1; jp <=npersons; jp++) {
		j = ncov + jp
		theta[j] = theta[j] * d[jp]
	}

	for (jf = 1; jf <=nfirms-1; jf++) {
		j = ncov + npersons + jf
		theta[j] = theta[j] * f[jf]
	}

	/** Writes estimation results in file cgout */

	fp = fopen("cgout_$STANDARDIZE_$TRANSFORMATION_$DEPENDENT","w");

	fputmatrix(fp, theta);

	fclose(fp);
}

void modcg(real vector b, real vector x, real scalar maxit, real scalar resid, 
		real matrix xd, 	real matrix xf, 
		real matrix df,
		real matrix dfp, 	real matrix dff,
		real scalar ncov, real scalar npersons,
		real scalar nfirms, real scalar ncells,
		real scalar ncoef, real scalar n) {

/* At the beginning x is the initial guess, x is then the approximate solution on output */

	real scalar eps
	real scalar info
	real scalar tol
	real scalar beta
	real scalar bnrm2
	real scalar itmax	

	real vector r , p , q

	printf("Starting Conjugate Gradient Algorithm\n");

	eps = 10e-15
	info = 0
	tol = resid
	beta = 0

	r = b
	p = J(ncoef, 1, 0) 

	bnrm2 = sqrt(b*b')

	if (bnrm2 == 0) {
		bnrm2 = 1
	}
	
	matvec (x, r, xd, xf, df, dfp, dff, ncov, npersons, nfirms, ncells, ncoef);

	r = b - r'

	itmax = maxit
	maxit = 0
	rnrm2 = sqrt(r*r')
	resid = rnrm2 / bnrm2

	printf("Iteration %f, norm of residual %f, relative error %f\n", maxit, rnrm2, resid);

	if (resid <= tol) return
	
	w = r

	rho = r*w'

	for (maxit = 1 ; maxit <= itmax ; maxit ++) {
		p = w' + beta * p
		matvec (p, q, xd, xf, df, dfp, dff, ncov, npersons, nfirms, ncells, ncoef)
		alpha = rho / (p'*q)
		x = x + alpha * p 
		r = r - alpha *q'
		rnrm2 = sqrt(r*r')
		resid = rnrm2 / bnrm2
		printf("Iteration %f, norm of residual %f, relative error %f\n", maxit, rnrm2, resid);
		if (rho < n*eps) return
		if (resid <= tol ) return
		rho1 = rho
		w = r
	 	rho = r*w'
		beta = rho/rho1
	}	

}

void matvec(real vector xin, 	real vector rout, 
		real matrix xd, 	real matrix xf,
		real matrix df,
		real matrix dfp, 	real matrix dff,
		real scalar ncov, real scalar npersons,
		real scalar nfirms, real scalar ncells,
		real scalar ncoef ) {

/* 	Computes the matrix vector product rout <- A*xin
	with A =(X D F)'(X D F)
	
	X'X in xxin(ncov,ncov), assumed identity from preconditioning
	X'D in xd(npersons, ncov)
	X'F in xf(nfirms, ncov)
	D'F in df(ncells), person and firm indices in dfp(ncells) and dff(ncells)
	D'D is identity from preconditioning dimension(npersons, npersons)
	F'F is identity from preconditioning dimension(nfirms, nfirms)

	The vectors X and R have three parts,
	The covariate effects 1:ncov
	The person effects ncov+1:ncov+npersons
	The firm effects ncov+npersons+1:ncov+npersons+nfirms 

*/

	real scalar i;

/* First the covariate effects */
	rout = xin
	rout[1 .. ncov] = cross(xd , (xin[ncov+1 .. ncov+npersons])) + rout[1 .. ncov]
	rout[1 .. ncov] 	 = cross(xf[1 .. nfirms - 1, 1 .. ncov ], (xin[ncov+npersons+1 .. ncov+npersons+nfirms-1]) ) + rout[1 .. ncov]

/* then the person effects */
	rout[ncov+1 .. ncov+npersons ] 			= xd * xin[1 .. ncov] + rout[ncov+1 .. ncov+npersons ] 

/* and finally the firm effects */
	rout[ncov+npersons+1 .. ncov+npersons+nfirms-1] 	= xf[1 .. nfirms - 1 , 1 .. ncov] * xin[1 .. ncov] + rout[ncov+npersons+1 .. ncov+npersons+nfirms-1]
	for ( i = 1 ; i <= ncells ; i++ ) {
		jperson 	= dfp[i] + ncov
		jfirm 	= dff[i] + ncov + npersons
	
		if (jfirm <= ncoef) {
			rout[jperson] 	= rout[jperson] 	+ xin[jfirm] * df[i]
			rout[jfirm] 	= rout[jfirm] 	+ xin[jperson] * df[i]
		}
	}

}



/* This function multiplies X's */
void xtprod(real vector persons, real vector firms, real matrix cov,
		real vector f, real vector d, 
		real vector s, real vector r,
		real scalar ncov, real scalar npersons,
		real scalar nfirms, real scalar ncells,
		real scalar n ) {

	/*
		Multiplies X's -> r
		X is of size ( ncov x n ) and s of size ( ncov )

	 	persons(n) 	: vector containing the $personid of the record 
		firms(n)	: vector containing the $firmid of the record
		cov(n,ncov)	: the covariates of the record
		s(n)		: the vector being multiplied
		f(nfirms) : vector containing  the number of times the firm
					appears in the data 
		d(npersons)	: vector containing the number of times the person
					appears in the data
		r(ncoef)	: the output vector
	*/

	real scalar i;
	real scalar tmp;
	real vector jp, jf;

	printf("XTPROD : Computing X'y\n");

	r = J(ncov + npersons + nfirms -1 ,1,0);

	for (i=1; i<=n; i++) {
		tmp = s[i]
		jp = persons[i]
		j = jp + ncov
		r[j] = r[j] + tmp * d[jp]
		jf = firms[i]
		if (jf != nfirms) {
			j = jf+ncov+npersons
			r[j] = r[j] + tmp * f[jf]

		}
	}		
	
	r[1..ncov] = cross(cov,s)

}

maincg(veccov);

end
