
program define a2reg
* A2REG uses mata, therefore requires version 9
version 9

* Syntax: a2reg depvar indepvars, individual(first fixed effect) unit(second fixed effect) save(output file)
syntax varlist(min=2), individual(varname) unit(varname) [prefix(name)] [save(name)]



* Checks whether all variables are numeric

foreach v in `varlist' {
	confirm numeric variable `v'
}

gettoken dependent varlist : varlist

tempvar indid unitid cellid

if ("`save'"=="") {
	tempfile save
}

if ("`prefix'"=="") {
	local prefix "effect"
}

* Creating sequenced ID variables

egen `indid'  = group(`individual')
egen `unitid' = group(`unit')
egen `cellid' = group(`indid' `unitid')

* Getting the size of the problem

sort `indid'
summarize `indid', meanonly 
local ninds = r(max) 

sort `unitid'
summarize `unitid', meanonly 
local nunits = r(max) 

sort `cellid'
summarize `cellid', meanonly 
local ncells = r(max) 

sort `indid' `unitid'

* Launches the Conjugate gradient estimation
* Note : the program fails when explanatory variables are not linearly independent

mata: maincg("`dependent'",tokens("`varlist'"),"`indid'","`unitid'", `ninds', `nunits', `ncells', "`save'")
mata: addresultscg("`dependent'",tokens("`varlist'"),"`indid'","`unitid'", `ninds', `nunits', `ncells', "`save'","`prefix'")

summarize indeffect, meanonly
local mean_ie = r(mean)

replace indeffect = indeffect - `mean_ie'  

end

mata

function maincg(string scalar dependent, string rowvector covariates,
				string scalar individualid, string scalar unitid,
				real scalar npupils, real scalar nschools,
				real scalar ncells, string scalar filerawoutput) {

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
	st_view(y=., 	., dependent)

	printf("Pupil IDs\n");
	st_view(pupil=., 	., individualid)

	printf("School IDs\n");
	st_view(school=., ., unitid)

	printf("Dependent variable : %s\n", dependent);



	printf("Problem Size ...\n");

	n 		    = 	length(y);
	ncov		=	cols(cov);
	ncoef		= 	npupils + nschools + ncov - 1;

	printf("%f observations, %f covariates, %f pupils, %f schools, %f cells\n",n,ncov,npupils,nschools,ncells);

	if ((length(y) != rows(cov)) | (length(y)!=length(pupil)) | (length(y)!=length(school)) ) {
		printf("Error : invalid matrix length\n");
		exit();
	}

	printf("Initializing normal form matrices\n");

	icell = 1;
	xd	= J(npupils,ncov,0);
	xf	= J(nschools,ncov,0);
	d	= J(npupils,1, 0);
	dfp	= J(ncells, 1, 0);
	dff	= J(ncells, 1, 0);
	df	= J(ncells, 1, 0);
	f	= J(nschools, 1, 0);

	printf("Starting loop\n");

	for ( i=1 ; i<=n ; i++ ) {

		if (i>1) {
			if (pupil[i] != pupil[i-1] || (school[i] != school[i-1] && pupil[i] == pupil[i-1])) {
				icell = icell +1;
			}
		}

		/* Construction de d */
		//printf("d ");
		jp = pupil[i]
		d[jp] = d[jp] + 1

		/* Construction de dfp */
		//printf("dfp ");
		dfp[icell] = jp

		/* Construction de f */
		//printf("f ");
		jf = school[i]
		f[jf] = f[jf] + 1

		/* Construction de dff */
		//printf("dff ");
		dff[icell] = jf

		/* Construction de df */
		//printf("df ");
		df[icell] = df[icell] + 1

		for ( j = 1; j <= ncov ; j++ ) {
		/* Construction de XD et XF */
	
			//printf("xd ");
			xd[jp,j] = xd[jp, j] + cov[i, j]
			//printf("xf ");
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

	for (i = 1; i<=npupils ; i ++) {
		d[i] = 1/sqrt(d[i]);
	}

	for (i = 1; i<=nschools ; i ++) {
		f[i] = 1/sqrt(f[i]);
	}

	for (j=1;j<=ncov;j++) {
		for (i = 1; i<=npupils ; i ++ ) {
			xd[i, j] = xd[i, j]*d[i];
		}
		for (i = 1; i<=nschools ; i ++ ) {
			xf[i, j] = xf[i, j]*f[i];
		}
	}

	for (i=1; i<=ncells; i++) {	
		df[i] = df[i]*d[dfp[i]]*f[dff[i]]
	}

	// Does (X D F)'y and puts it in b 
	xtprod (pupil, school, cov[1 .. n,1..ncov], f, d, y, b, ncov, npupils, nschools, ncells, n);
	b = b'

	resid = 1.0e-7
	maxit = 1000

	printf("Scaling theta for preconditioning\n");

	u
	theta[1 .. ncov]
	
	theta[1..ncov] = u * theta[1..ncov]

	printf("Pupil effects\n");

	for (jp = 1; jp<=npupils; jp++) {
		j = ncov + jp
		theta[j]=theta[j]/d[jp]	
	}

	printf("School effects\n");

	for (jf = 1; jf<=nschools-1; jf++) {
		j = ncov + npupils +jf
		theta[j] = theta[j]/f[jf]
	}
	
	printf("Beginning Conjugate Gradient Iterations\n");

	modcg(b,theta,maxit,resid, xd, xf, df, dfp, dff, ncov, npupils, nschools, ncells, ncoef, n);

	printf("Transforming theta back to original scale\n");

	printf("Size of u (%f x %f) and of theta (%f x %f)\n", rows(u), cols(u), rows(theta), cols(theta));

	theta[1 .. ncov] = luinv(u) * theta[1 .. ncov]

	for (jp = 1; jp <=npupils; jp++) {
		j = ncov + jp
		theta[j] = theta[j] * d[jp]
	}

	for (jf = 1; jf <=nschools-1; jf++) {
		j = ncov + npupils + jf
		theta[j] = theta[j] * f[jf]
	}

	/** Writes estimation results in file cgout */

	fp = fopen(filerawoutput,"w");

	fputmatrix(fp, theta);

	fclose(fp);
}

void modcg(real vector b, real vector x, real scalar maxit, real scalar resid, 
		real matrix xd, 	real matrix xf, 
		real matrix df,
		real matrix dfp, 	real matrix dff,
		real scalar ncov, real scalar npupils,
		real scalar nschools, real scalar ncells,
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
	
	matvec (x, r, xd, xf, df, dfp, dff, ncov, npupils, nschools, ncells, ncoef);

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
		matvec (p, q, xd, xf, df, dfp, dff, ncov, npupils, nschools, ncells, ncoef)
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
		real scalar ncov, real scalar npupils,
		real scalar nschools, real scalar ncells,
		real scalar ncoef ) {

/* 	Computes the matrix vector product rout <- A*xin
	with A =(X D F)'(X D F)
	
	X'X in xxin(ncov,ncov), assumed identity from preconditioning
	X'D in xd(npupils, ncov)
	X'F in xf(nschools, ncov)
	D'F in df(ncells), person and firm indices in dfp(ncells) and dff(ncells)
	D'D is identity from preconditioning dimension(npupils, npupils)
	F'F is identity from preconditioning dimension(nschools, nschools)

	The vectors X and R have three parts,
	The covariate effects 1:ncov
	The pupil effects ncov+1:ncov+npupils
	The school effects ncov+npupils+1:ncov+npupils+nschools 

*/

	real scalar i;

/* First the covariate effects */
	rout = xin
	rout[1 .. ncov] = cross(xd , (xin[ncov+1 .. ncov+npupils])) + rout[1 .. ncov]
	rout[1 .. ncov] 	 = cross(xf[1 .. nschools - 1, 1 .. ncov ], (xin[ncov+npupils+1 .. ncov+npupils+nschools-1]) ) + rout[1 .. ncov]

/* then the pupil effects */
	rout[ncov+1 .. ncov+npupils ] 			= xd * xin[1 .. ncov] + rout[ncov+1 .. ncov+npupils ] 

/* and finally the school effects */
	rout[ncov+npupils+1 .. ncov+npupils+nschools-1] 	= xf[1 .. nschools - 1 , 1 .. ncov] * xin[1 .. ncov] + rout[ncov+npupils+1 .. ncov+npupils+nschools-1]
	for ( i = 1 ; i <= ncells ; i++ ) {
		jpupil 	= dfp[i] + ncov
		jschool 	= dff[i] + ncov + npupils
	
		if (jschool <= ncoef) {
			rout[jpupil] 	= rout[jpupil] 	+ xin[jschool] * df[i]
			rout[jschool] 	= rout[jschool] 	+ xin[jpupil] * df[i]
		}
	}

}



/* This function multiplies X's */
void xtprod(real vector pupils, real vector schools, real matrix cov,
		real vector f, real vector d, 
		real vector s, real vector r,
		real scalar ncov, real scalar npupils,
		real scalar nschools, real scalar ncells,
		real scalar n ) {

	/*
		Multiplies X's -> r
		X is of size ( ncov x n ) and s of size ( ncov )

	 	pupils(n) 	: vector containing the pupilid of the record 
		schools(n)	: vector containing the schoolid of the record
		cov(n,ncov)	: the covariates of the record
		s(n)		: the vector being multiplied
		f(nschools) : vector containing  the number of times the school
					appears in the data 
		d(npupils)	: vector containing the number of times the pupil
					appears in the data
		r(ncoef)	: the output vector
	*/

	real scalar i;
	real scalar tmp;
	real vector jp, jf;

	printf("XTPROD : Computing X'y\n");

	r = J(ncov + npupils + nschools -1 ,1,0);

	for (i=1; i<=n; i++) {
		tmp = s[i]
		jp = pupils[i]
		j = jp + ncov
		r[j] = r[j] + tmp * d[jp]
		jf = schools[i]
		if (jf != nschools) {
			j = jf+ncov+npupils
			r[j] = r[j] + tmp * f[jf]

		}
	}		
	
	r[1..ncov] = cross(cov,s)

}

function addresultscg(string scalar dependent, string rowvector covariates,
				string scalar individualid, string scalar unitid,
				real scalar npupils, real scalar nschools,
				real scalar ncells, string scalar filerawoutput,

				string scalar prefix) {




	real scalar ncov;



	ncov = length(covariates);	



	fp = fopen(filerawoutput,"r")

	params = fgetmatrix(fp);

	

	betas = params[1 .. ncov];

	pupileffects = params[ncov+1 .. ncov + npupils];

	schooleffects = params[ncov+npupils+1 .. ncov + npupils + nschools -1];



	st_addvar("double", "indeffect");

	st_addvar("double", "uniteffect");	

	st_view(data, . ,(individualid, unitid, "indeffect","uniteffect"));



	n = rows(data);

	printf(" %f observations ", n);



	for (i = 1; i<= n ; i ++) {

		data[i,3] = pupileffects[data[i,1]];

		if (data[i,2] != nschools) {

			data[i,4] = schooleffects[data[i,2]];

		} else {

			data[i,4] = 0;

		}

	}

	

	for (i = 1; i<= ncov; i++) {

		stata(sprintf("gen %s_%s = %f", prefix, covariates[i], betas[i]));

		stata(sprintf("label variable %s_%s %s Coefficient for %s %s ", prefix, covariates[i], char(34),covariates[i],char(34)));

	}

}



end


