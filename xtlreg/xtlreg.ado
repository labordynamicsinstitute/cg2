/*
 *
 * XTLREG: Estimates models of two-way fixed effects with one lagged fixed effect
 *
 */
 





#delimit ;


cap program drop xtlreg;

/*
 * Syntax: 
 * 	xtlreg y x1 x2 , individual(personid) unit(unitid) lambda(lambda) time(timeid) sorted
 *
 *	xtlreg estimates a model with y as the dependent variable, x1, x2 as the explanatory variables
 *	an individual effect personid, a unit effect unitid, lagged in the second period,
 *	and with a discounting factor lambda, time period indexed by timeid
 *
 */

 
 
/*
 *
 * eclass: returns estimation results
 * sortpreserve: preserves the sorting order
 *
 */

program define xtlreg, eclass sortpreserve;
/* XTLREG uses mata, therefore requires version 9 */
version 9;

global 	AUTHOR	=	"Amine Ouazad";
global	VERSION	=	"0.1";

syntax varlist(min=2), individual(varname) unit(varname) 
	lambda(real) time(varname) sorted
	[indeffect(name)] [uniteffect(name)]
	[pastuniteffect(name)]
	[resid(name)] [xb(name)] [largestgroup];

display "xtlreg program version $VERSION - Estimation of two-way fixed effects with past unit effects";

/* Checks whether all variables are numeric */
display "Checks whether all variables are numeric";

foreach v in `varlist' {;
	confirm numeric variable `v';
};

if ("`largestgroup'"!="") {;
	quietly{;
	tempvar group group_count;
	a2group, individual(`individual') unit(`unit') groupvar(`group');
	egen `group_count' = count(`group'), by(`group');
	summarize `group_count';
	local max = r(max);
	summarize `group' if `group_count' == `max';
	local largest_group = r(min);
	drop if `group' != `largest_group';
	};
};

if ("`indeffect'"!="") {;
	confirm new variable `indeffect';
};

if ("`indeffect'"=="") {;
	tempvar indeffect;
};

if ("`uniteffect'"!="") {;
	confirm new variable `uniteffect';
};

if ("`uniteffect'"=="") {;
	tempvar uniteffect;
};

if ("`pastuniteffect'"!="") {;
	confirm new variable `pastuniteffect';
};

if ("`pastuniteffect'"=="") {;
	tempvar pastuniteffect;
};

if ("`resid'"!="") {;
	confirm new variable `resid';
};

if ("`resid'"=="") {;
	tempvar resid;
};

if ("`xb'"!="") {;
	confirm new variable `xb';
};

if ("`xb'"=="") {;
	tempvar xb;
};

marksample touse, novarlist strok;
keep if `touse';

gettoken dependent varlist : varlist;

tempvar indid unitid pastunitid timeid cellid_uup cellid_iu cellid_iup;
tempfile save;

/*
 *
 * Creating sequenced ID variables
 *
 */

 display "Creating sequenced ID variables";

egen `indid'  = group(`individual');
egen `unitid' = group(`unit');
egen `timeid' = group(`time');

tempfile temp_ids temp_data;
save "`temp_data'", replace;

keep `indid' `timeid' `unitid';
duplicates drop `indid' `timeid', force;
tsset `indid' `timeid';
gen `pastunitid' =l.`unitid';
sort `indid' `timeid';
save "`temp_ids'", replace;

use "`temp_data'";
sort `indid' `timeid';
merge `indid' `timeid' using "`temp_ids'", keep(`indid' `timeid' `pastunitid');
drop _merge;

egen `cellid_iu' = group(`indid' `unitid');

sort `unitid' `pastunitid';
egen `cellid_uup' = group( `unitid' `pastunitid' ) ;

sort `indid' `pastunitid';
egen `cellid_iup' = group( `indid'  `pastunitid' );

summarize `indid', meanonly ;
local nind = r(max) ;

summarize `unitid', meanonly ;
local nunit = r(max) ;

summarize `cellid_iu', meanonly ;
local ncells_iu = r(max) ;

summarize `cellid_uup', meanonly ;
local ncells_uup = r(max) ;

summarize `cellid_iup', meanonly ;
local ncells_iup = r(max) ;

sort `indid' `unitid';

/* Launches the Conjugate gradient estimation */
/* Note : the program fails when explanatory variables are not linearly independent */

display " Dependent variable 	: `dependent' ";
display " Covariates		: `varlist' ";
display " Individual ID 	: `individual' ";
display " Unit ID 		: `unit' ";
display " Time 			: `time' ";

mata: maincgxtl("`dependent'",tokens("`varlist'"),"`indid'","`unitid'", "`pastunitid'", `lambda', `nind', `nunit', `ncells_iu', `ncells_uup', `ncells_iup',"`save'");

quietly {;

gen double `indeffect' = .;
gen double `uniteffect' = .;
gen double `pastuniteffect' = .;
gen double `xb' = .;

mata: addresultscgxtl("`dependent'",tokens("`varlist'"),"`indid'","`unitid'", "`pastunitid'", `nind', `nunit', "`save'","`indeffect'","`uniteffect'","`pastuniteffect'","`xb'");

gen `resid' = `dependent' - `xb' - `indeffect' - `uniteffect' - `lambda' * `pastuniteffect' if `pastuniteffect' != .;
replace `resid' = `dependent' - `xb' - `indeffect' - `uniteffect'  if `pastuniteffect' == .;

summarize `indeffect', meanonly;
scalar mean_ie = r(mean);
scalar stdev_ie = r(sd);

summarize `uniteffect', meanonly;
scalar mean_ue = r(mean);
scalar stdev_ue = r(sd);

summarize `pastuniteffect', meanonly;
scalar mean_upe = r(mean);
scalar stdev_upe = r(sd);

summarize `xb', meanonly;
scalar mean_xb = r(mean);
scalar stdev_xb = r(sd);

summarize `resid', meanonly;
scalar mean_resid = r(mean);
scalar stdev_resid = r(sd);

scalar RSS = mean_resid * _N;

replace `indeffect' = `indeffect' - mean_ie;

matrix means = (mean_ie, 	mean_ue, 	mean_upe, 	mean_xb, 	mean_resid );
matrix colnames means = Individual_Effect Unit_Effect Past_Unit_Effect XB Residual;

matrix sdevs = (stdev_ie, 	stdev_ue, 	stdev_upe, 	stdev_xb, 	stdev_resid);
matrix colnames sdevs = Individual_Effect Unit_Effect Past_Unit_Effect XB Residual;

matrix accum R = y indeffect uniteffect resid, nocons dev;
matrix R = corr(R);
matrix colnames R = Dependent Individual_Effect Unit_Effect Residual;
matrix rownames R = Dependent Individual_Effect Unit_Effect Residual;

ereturn scalar nobs 	= _N;
ereturn scalar nind     = `nind';
ereturn scalar nunit	= `nunit';
ereturn scalar constant	= mean_ie;

ereturn matrix beta	betas;
ereturn matrix M		means;
ereturn matrix SD		sdevs;
ereturn matrix corrtable R;
};
end;

#delimit cr



mata

function maincgxtl(string scalar dependent, string rowvector covariates,
			string scalar individualid, string scalar unitid, string scalar pastunitid,
			real scalar lambda,
			real scalar nind, real scalar nunit, 
			real scalar ncells_iu, real scalar ncells_uum, real scalar ncells_ium,
			string scalar filerawoutput) {

/*

 
	STEPS OF THE CONJUGATE GRADIENT ESTIMATION
    ==========================================
    
    1. Build matrices or representations of the matrices :
		X'X X'D X'F X'F-1 (brutal multiplication by cross)
		D'D F'F F-1'F-1' (three loops)
		D'F D'F-1 F'F-1  (building cells - sparse representation)
		and X'Y (xtprod)

	2. Precondition:
		Perform the Cholesky decomposition of X'X = U'U
		Transform matrices off the diagonal

	3. Start conjugate gradient iterations
		This requires knowing how to multiply by A (modcg)

*/

/*************** DATA DECLARATION */

	real scalar n, ncov, ncoef;
	
	real scalar i, j, jp, jf, jfm; 		// Used for storing person and firm indices

	real vector y;				// The vector of observations (Y)					-	Size : n
	real matrix cov;			// The covariate matrix X (view)					- 	Size : n*ncov
	real vector individual;			// The vector of individual indices i				- 	Size : n
	real vector unit;			// The vector of unit indices J(i,t)				- 	Size : n
	real vector pastunit;			// The vector of past unit indices J(i,t-1)			-	Size : n

	real vector d,   delta;          	// Sparse representation of matrix D'D F'F F'F-1		-	Size : nind + nunit-1

	real vector df,  dfp , dff;        	// Sparse representation of matrix D'F				-	Size : 3*ncell_iu
	real vector dfm, dfmp, dfmfm;      	// Sparse representation of matrix D'F-1				- 	Size : 3*ncell_uum
	real vector ffm, ffmf, ffmfm;      	// Sparse representation of matrix F'F-1				- 	Size : 3*ncell_ium
		
	real matrix xx, xd, xfl;	     	// Matrices X'X X'D X'phi_lambda 					-	Size : ncov²+ncov*nind+ncov*(nunit-1)

	real matrix u ;			// Cholesky decomposition of X'X = U'U				- 	Size : ncov²/2
	real matrix upinv;			// inverse of U'

	real vector b ; 			// Right hand side of the normal form equation of the OLS	-	Size : n

	real vector theta ; 			// The parameter vector 						- 	Size : ncoef = ncov + nind + nunit-1 

	real scalar icell_iu , icell_uum, icell_ium; 	// Index used to scan through the cells of the dataset

	real scalar maxit; 			// Maximum number of conjugate gradient iterations
	real scalar resid;			// Residual

	real scalar fp;				// File descriptor to save the parameter vector

	real scalar fmissing;			// Value of missing firm	
	
	

/***************** GETTING VIEWS ON THE DATA */

	st_view(cov=., 		., covariates);
	st_view(y=., 		., dependent);
	st_view(individual=., 	., individualid);
	st_view(unit=., 	., unitid);
	st_view(pastunit=., 	., pastunitid);

	fmissing = missingof(pastunit);

	n     = length(y);
	ncov  = cols(cov);
	ncoef = ncov + nind + nunit -1  ;

	printf("%f observations, %f covariates, %f individuals, %f units\n",n,ncov,nind,nunit);
	printf("%f D'F cells, %f F'F-1 cells, %f D'F-1 cells \n", ncells_iu, ncells_uum, ncells_ium);
	printf("A total of %f coefficients\n",	ncoef);	
	printf("Lambda is : %f\n", lambda);

	theta = J( ncoef, 1, 0 );

	xd	= J( ncov,  	nind, 	 0);
	xfl	= J( ncov, 	nunit-1, 0);

	d     = J( nind,        1,    0);
	delta = J( nunit-1,	1,    0);

	dfp   = J( ncells_iu,   1,    0);
	dff   = J( ncells_iu,   1,    0);
	df    = J( ncells_iu,   1,    0);

	ffmf  = J( ncells_uum,   1,    0);
	ffmfm = J( ncells_uum,   1,    0);
	ffm   = J( ncells_uum,   1,    0);

	dfmp  = J( ncells_ium,   1,    0);
	dfmfm = J( ncells_ium,   1,    0);
	dfm   = J( ncells_ium,   1,    0);

	b	= J( ncoef,		 1,	 0);



	if ((length(y) != rows(cov)) | (length(y)!=length(individual)) | (length(y)!=length(unit)) | (length(y)!=length(pastunit)) ) {
		printf("Error : invalid matrix length\n");
		exit();
	}


	printf("Building X'X... ");
	xx = cross(cov,cov)
	xx

	printf("... size %g x %g\n", rows(xx), cols(xx));

	printf("Building D'D, X'D, X'Phi_lambda\n");
	for ( i=1; i<=n ; i++ ) {
		//printf("* \t");
		//printf("i : %g \t",i);
		jp 		= individual[i] ;
		//printf("jp : %g \n",jp);
		d[jp] 	= d[jp] + 1 ;
		
		jf 		= unit[i] ;
		//printf("jf : %g \n",jf);
		if (jf != nunit) {
			assert(jf < nunit);
			delta[jf] = delta[jf] + 1 ;
		}

		jfm 		= pastunit[i] ;
		//printf("jfm : %g \n", jfm);

		if (jfm != fmissing && jfm >= 1 && jfm != nunit) {
			delta[jfm] = delta[jfm] + lambda 		* lambda ;
			assert(jfm < nunit);
		}

		for ( j = 1; j <= ncov ; j++ ) {
			/* Construction de XD XF + lambda*XF-1 				*/
			/* Cela suppose que xd xf et xfm soient initialisés à 0 	*/
			xd[j,jp]   = xd[j,jp] 	+ cov[i, j] ;
			if (jf  != nunit) {
				xfl[j,jf]  = xfl[j,jf ] + 		cov[i, j];
			}
			if (jfm != fmissing && jfm >= 1 && jfm != nunit) {
				xfl[j,jfm] = xfl[j,jfm] + lambda  *	cov[i, j];
			}
		}

	}

/************ BUILDING SPARSE REPRESENTATIONS OF D'F F'F-1 D'F-1 */

// D'F

	printf("Building D'F\n");
	stata(sprintf("sort %s %s",individualid,unitid));

	icell_iu = 1;
	for ( i=1 ; i<=n ; i++ ) {
		if ( i > 1 ) {
			if (individual[i] != individual[i-1] || (unit[i] != unit[i-1] && individual[i] == individual[i-1])) {
				icell_iu = icell_iu +1;
			}
		}
		dfp[icell_iu] = individual[i];
		dff[icell_iu] = unit[i];
		df[icell_iu]  = df[icell_iu] + 1;

	}

	assert(icell_iu == ncells_iu);

// D'F-1

	printf("Building D'F-1\n");	
	stata(sprintf("sort %s %s",individualid,pastunitid));
	
	icell_ium = 1;
	real scalar notfirst ;
	notfirst =0;
	for ( i = 1  ; i <= n ; i++ ) {
		if (pastunit[i] != fmissing) {
			if ( i > 1 && notfirst) {
				if ((individual[i] != individual[i-1]) || (pastunit[i] != pastunit[i-1] && individual[i] == individual[i-1]))   {
					icell_ium = icell_ium + 1;
				}
			}
			notfirst = 1;
			if (icell_ium > ncells_ium) {
				printf("Too many cells\n");
				printf("i : %g \t individual[i] : %g \t individual[i-1] : %g \t pastunit[i] : %g \t icell_ium : %g \n", i, individual[i], individual[i-1], pastunit[i], icell_ium);
				assert(0);
			}
			dfmp[icell_ium]  = individual[i];
			dfmfm[icell_ium] = pastunit[i];
			dfm[icell_ium]   = dfm[icell_ium] + 1;		
		}
	}
	assert(icell_ium == ncells_ium);

// F'F-1

	printf("Building F'F-1\n");	
	stata(sprintf("sort %s %s",unitid,pastunitid));
	
	icell_uum = 1;
	notfirst = 0;
	for ( i=1 ; i<=n ; i++ ) {
		if (pastunit[i] != fmissing) {
			if ( i > 1 && notfirst ) {
				if ((unit[i] != unit[i-1]) || (pastunit[i] != pastunit[i-1] && unit[i] == unit[i-1])) {
					icell_uum = icell_uum + 1;
				}
			}
			notfirst = 1;
			assert(icell_uum <= ncells_uum);
			ffmf[icell_uum]  = unit[i];
			ffmfm[icell_uum] = pastunit[i];
			ffm[icell_uum]   = ffm[icell_uum] + 1;
			
			/* Adds elements to the diagonal of Phi_lambda' Phi_lambda */
			if (ffmf[icell_uum] == ffmfm[icell_uum] && ffmfm[icell_uum]!=nunit) {
				delta[ffmf[icell_uum]] = delta[ffmf[icell_uum]] + 2.0 * lambda;
				/* Notice that Phi_Lambda'Phi_Lambda = F'F + lambda^2 F_-1 ' F_-1 
					+ lambda * F_-1'F + lambda * F'F_-1
					so that actually a diagonal element adds 2 times lambda */
			}
		}
	}

	assert(icell_uum == ncells_uum);


/********** PRECONDITIONING */

	/*
	 *
	 * Preconditioning is done by multiplying the diagonal so that the diagonal is filled
	 * with ones
	 *
	 */

	printf("Starting preconditioning\n");

	u = cholesky(xx)';
	upinv = luinv(u');

	cov	= cov * upinv' ;

	xd 	= upinv * xd;
	xfl 	= upinv * xfl;

	for (i = 1; i<=nind ; i ++) {
		d[i] = 1/sqrt(d[i]);
	}

	for (i = 1; i<=nunit-1 ; i ++) {
		delta[i] = 1/sqrt(delta[i]);
	}

	for ( i = 1 ; i <= ncov ; i++) {
		for (j = 1; j <= nind ; j ++ ) {
			xd[i, j]  = xd[i, j]	*	d[j];
		}
		for (j = 1; j <= nunit-1 ; j ++ ) {
			xfl[i, j] = xfl[i, j]	*	delta[j];
		}
	}


	for (i=1; i<=ncells_iu; i++) {
		if (dff[i] < nunit) {	
			df[i] = df[i]*d[dfp[i]]*delta[dff[i]]
		}
	}

	for (i=1; i<=ncells_ium; i++) {
		if (dfmfm[i] < nunit) {
			dfm[i] = dfm[i]*d[dfmp[i]]*delta[dfmfm[i]]
		}
	}

	for (i=1; i<=ncells_uum; i++) {	
		if (ffmf[i] < nunit && ffmfm[i] < nunit ) {
			ffm[i] = ffm[i]*delta[ffmf[i]]*delta[ffmfm[i]]
		}
	}

	xtprodxtl (individual, unit, pastunit, cov, y, d, delta, b, lambda,
			n, ncov, nind, nunit );

	resid = 1.0e-7
	maxit = 1000


	printf("Beginning Conjugate Gradient Iterations\n");



	modcgxtl(b, theta, maxit, resid, xd, xfl,
		df, dfp, dff, ncells_iu,
		dfm, dfmp, dfmfm,	ncells_ium,
		ffm, ffmf, ffmfm,	ncells_uum,
		lambda,
		ncov, nind, nunit, ncoef, n);

	printf("Transforming theta back to original scale\n");


	theta[1 .. ncov] = luinv(u) * theta[1 .. ncov];

	for (jp = 1; jp <= nind; jp++) {
		j = ncov + jp
		theta[j] = theta[j] * d[jp]
	}


	for (jf = 1; jf < nunit; jf++) {
		j = ncov + nind + jf
		theta[j] = theta[j] * delta[jf]
	}

	//printf("Saving to file %s\n",filerawoutput);
	fp = fopen(filerawoutput,"w");

	fputmatrix(fp, theta);

	fclose(fp);


}


void xtprodxtl(real vector individual, real vector unit, real vector pastunit, real matrix cov, // the X = (beta theta psi psi-1) matrix
		real vector y, real vector d, real vector delta, real vector r,
		real scalar lambda,
		real scalar n, real scalar ncov, real scalar nind, real scalar nunit
		)

///////////////////////////////////////////////////////////
//
// MULTIPLIES X's = r
//
///////////////////////////////////////////////////////////

{

/////////////////////////////////
//		( X's  ) 		 //
// X's  = 	( D's ) 		 //
//		( Phi_Lambda's   ) //
/////////////////////////////////

	real scalar i, jind, junit, jpastunit, tmp, j, k, l;	
	real scalar fmissing;
	real scalar ncoef;

	fmissing = missingof(pastunit);

	ncoef = ncov + nind + nunit - 1 ;
	r = J( ncoef , 1 , 0);

// First : beta's

	r[1..ncov]  = cross(cov,y)

// Then the theta's psi's and psi-1's

	for ( i = 1 ; i <= n ; i ++ ) {
		
		tmp 		= y[i];
		jind 		= individual[i];
		junit 	= unit[i];
		jpastunit 	= pastunit[i];

		j = ncov + jind
		k = ncov + nind + junit
		l = ncov + nind + jpastunit 
		
		r[ j ] = r[ j ] +	d[jind] 	  * tmp;

		if ( junit <= nunit - 1) {
			r[ k ] = r[ k ] +	delta[junit] 	  * tmp;    
		}

		if ( jpastunit != fmissing && jpastunit <=nunit - 1) {   // Only if the past unit id is less than the maximum and not missing
			r[ l ]  = r[ l ] + lambda * delta[jpastunit] * tmp;
		}

	}

}

void modcgxtl(	real vector b, real vector x, real scalar maxit, real scalar resid, 
		real matrix xd, real matrix xfl, 
		real matrix df, real matrix dfp, 	real matrix dff,	real scalar ncells_iu,
		real matrix dfm, real matrix dfmp, real matrix dfmfm,	real scalar ncells_ium,
		real matrix ffm, real matrix ffmf, real matrix ffmfm,	real scalar ncells_uum,
		real scalar lambda,
		real scalar ncov, real scalar nind, real scalar nunit,
		real scalar ncoef, real scalar n) {

/* At the beginning x is the initial guess, x is then the approximate solution on output */

	real scalar eps
	real scalar tol
	real scalar beta
	real scalar bnrm2
	real scalar itmax	
	real scalar rnrm2, w, rho, alpha, rho1

	real vector r , p , q

	eps = 10e-15
	tol = resid
	beta = 0

	r = b
	p = J(ncoef, 1, 0) 

	bnrm2 = sqrt(b'*b)

	if (bnrm2 == 0) {
		bnrm2 = 1
	}
	
	matvecxtl (	x, r, 
			xd, 	xfl,	
			df, 	dfp,	dff, 		ncells_iu,
			dfm,	dfmp,	dfmfm,	ncells_ium,
			ffm,  ffmf,	ffmfm,	ncells_uum,
			lambda,
			ncov, nind, nunit, ncoef );

	r = b - r

	itmax = maxit
	maxit = 0
	rnrm2 = sqrt(r'*r)
	resid = rnrm2 / bnrm2

	printf("Iteration %f, norm of residual %f, relative error %f\n", maxit, rnrm2, resid);

	if (resid <= tol) return
	
	w = r

	rho = r'*w

	for (maxit = 1 ; maxit <= itmax ; maxit ++) {
		p = w + beta * p
		matvecxtl (p, q, xd, 	xfl,	
			df, 	dfp,	dff, 		ncells_iu,
			dfm,	dfmp,	dfmfm,	ncells_ium,
			ffm,  ffmf,	ffmfm,	ncells_uum,
			lambda,
			ncov, nind, nunit, 	ncoef);
		alpha = rho / (p'*q)
		x = x + alpha * p 
		r = r - alpha * q
		rnrm2 = sqrt(r'*r)
		resid = rnrm2 / bnrm2
		printf("Iteration %f, norm of residual %f, relative error %f\n", maxit, rnrm2, resid);
		if (rho < n*eps) return
		if (resid <= tol ) return
		rho1 = rho
		w = r
	 	rho = r'*w
		beta = rho/rho1
	}	

}


void matvecxtl (
		real vector xin, real vector rout,
		real matrix xd,   real matrix xfl,						
		real matrix df,   real matrix dfp,	real matrix dff,		real scalar ncells_iu,
		real matrix dfm,  real matrix dfmp,	real matrix dfmfm,	real scalar ncells_ium,
		real matrix ffm,  real matrix ffmf,	real matrix ffmfm,	real scalar ncells_uum,
		real scalar lambda,
		real scalar ncov, 
		real scalar nind, 
		real scalar nunit,
		real scalar ncoef
		)

////////////////////////////////////////////////////////
//                                                 	//
// MATVEC computes rout = A*xin                    	//
//                                                 	//
// with A = ( X D Phi_lambda )' ( X D F Phi_Lambda ) 	//
//								     	//
////////////////////////////////////////////////////////

{
	
	//printf("Matvecxtl\n");

	real scalar i, jind, junit, jpastunit;

////////////////////////////////////////////////////////////////////
//		( beta  		+ X'D theta  	+ X'Phil psi   )	//
//  AX =    ( theta 		+ D'X beta   	+ D'Phil psi   )	//
//		( Phil'Phil psi   + Phil'X beta   	+ Phil'D theta )	//
// 										     	//
//	X'D, X'F, X'F-1 are stored as full matrices.               	//
//    Other matrices are stored as sparse matrices.		     	//
////////////////////////////////////////////////////////////////////

	rout = xin;

// First line
	//printf("First line\t");
	rout[1 .. ncov] = xd*xin[ncov+1 .. ncov+nind] 								+ rout[1 .. ncov]
	rout[1 .. ncov] = xfl*xin[ncov+nind+1 .. ncov+nind+nunit-1] 					+ rout[1 .. ncov]

// Second line without sparse matrices : adds X'D beta
	//printf("Second line\t");
	rout[ncov+1 .. ncov+nind ] = cross(xd, xin[1 .. ncov]) 						+ rout[ncov+1 .. ncov+nind] 

// Third line without sparse matrices  : puts F'X beta
	//printf("Third line\t");
	rout[ncov+nind+1 .. ncov+nind+nunit-1] = rout[ncov + nind + 1 .. ncov + nind + nunit-1] 	+ cross(xfl, xin[1 .. ncov]) 	

// Adding the contents of D'F
	//printf("Adds the contents of D'F\t");
	for (i = 1; i<= ncells_iu; i++) {
		jind  = ncov + dfp[i];
		junit = ncov + nind + dff[i];

		if (dff[i] != nunit) {
			rout[jind]  = xin[junit]*df[i] + rout[jind]
			rout[junit] = xin[jind] *df[i] + rout[junit]
		}
	}

// Adding the contents of D'F-1
	//printf("Adds the contents of D'F-1\t");
	for (i = 1 ; i <= ncells_ium ; i++ ) {
		jind 		= ncov + dfmp[i];
		jpastunit 	= ncov + nind + dfmfm[i];
		
		if (dfmfm[i] != nunit) {
			rout[jind] 		= rout[jind] 	+ lambda * xin[jpastunit]	* dfm[i]
			rout[jpastunit] 	= rout[jpastunit] + lambda * xin[jind] 		* dfm[i]
		}
	} 
	
// Adding the contents of F'F-1
	//printf("Adds the contents of F'F-1\n");
	for (i = 1 ; i <= ncells_uum ; i++ ) {
		junit		= ncov + nind + ffmf[i] ;
		jpastunit	= ncov + nind + ffmfm[i] ;

		if (ffmf[i] != nunit && ffmfm[i] != nunit && junit!=jpastunit) {
			rout[junit]		= rout[junit] 	+ lambda * xin[jpastunit] * ffm[i]
			rout[jpastunit]	= rout[jpastunit] + lambda * xin[junit] 	  * ffm[i]
		}
	}

}

function addresultscgxtl(string scalar dependent, string rowvector covariates,
				string scalar individualid, string scalar unitid, string scalar pastunitid,
				real scalar nind, real scalar nunit, 
				string scalar filerawoutput,
				string scalar indeffect_name,
				string scalar uniteffect_name,
				string scalar pastuniteffect_name,
				string scalar xb_name) {


	real scalar ncov;
	real scalar fp;
	real vector betas, indeffect, uniteffect;
	real vector params;
	real matrix data;
	real scalar n,i;
	string scalar indeffect_name, uniteffect_name, pastuniteffect_name;
	string scalar cmd_listcovnames;
	string scalar cmd_listcovvalues;

	ncov = length(covariates);
	
	fp = fopen(filerawoutput,"r");
	params = fgetmatrix(fp);
	fclose(fp);
	
	betas 		= params[1 		      	.. ncov				];
	indeffect 		= params[ncov + 1 		.. ncov + nind			];
	uniteffect 		= params[ncov + nind + 1 	.. ncov + nind + nunit -1	];

	st_view(data, . ,(individualid, unitid, pastunitid, indeffect_name, uniteffect_name, pastuniteffect_name));

	n = rows(data);

	for (i = 1; i<= n ; i ++) {
		data[i,4] = indeffect[data[i,1]];
		if (data[i,2] != nunit) {
			data[i,5] = uniteffect[data[i,2]];
		} else {
			data[i,5] = 0;
		}
		if (data[i,3] != nunit && data[i,3] != .) {
			data[i,6] = uniteffect[data[i,3]];
		} else if (data[i,3] == nunit) {
			data[i,6] = 0;
		} else {
			data[i,6] = .;
		}
	}

	cmd_listcovnames	=	"";
	cmd_listcovvalues	= 	"";

	for (i = 1; i<= ncov; i++) {
		cmd_listcovnames 	= sprintf("%s %s",cmd_listcovnames, covariates[i]);
		cmd_listcovvalues	= sprintf("%s %g,",cmd_listcovvalues, betas[i]);
		stata(sprintf("quietly replace %s = %s + %g * %s",xb_name,xb_name,betas[i], covariates[i]));
	}

	stata(sprintf("matrix input betas = (%s)",cmd_listcovvalues));
	stata(sprintf("matrix colnames betas = %s",cmd_listcovnames));
	
}

end

