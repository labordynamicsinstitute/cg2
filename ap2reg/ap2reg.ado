program define ap2reg, eclass
* A2REG uses mata, therefore requires version 9
version 9

global 	AUTHOR	=	"Amine Ouazad"
global	VERSION	=	"0.1"

* Syntax: ap2reg depvar indepvars, individual(first fixed effect) unit(second fixed effect) time(time variable) save(output file)
syntax varlist(min=2), individual(varname) unit(varname) time(varname) [prefix(name)] [save(name)]


display "ap2reg program version $VERSION - Estimation of two-way fixed effects with past unit effects"

quietly {

cap drop indeffect uniteffect pastuniteffect xb resid

* Checks whether all variables are numeric

foreach v in `varlist' {
	confirm numeric variable `v'
}

gettoken dependent varlist : varlist

tempvar timeid indid unitid pastunitid cellid_iu cellid_uup cellid_iup t

if ("`save'"=="") {
	tempfile save
}

if ("`prefix'"=="") {
	local prefix "effect"
}

* Creating sequenced ID variables

egen `timeid' = group(`time')
egen `indid'  = group(`individual')
egen `unitid' = group(`unit')

egen `cellid_iu' = group(`indid' `unitid')

sort `indid' `timeid'
gen `t' = .
replace `t' = `unitid'[_n - 1] if `indid'[_n-1] == `indid'[_n] 
egen `pastunitid' = group(`t')

sort `unitid' `pastunitid'
egen `cellid_uup' = group( `unitid' `pastunitid' ) 

sort `indid' `pastunitid'
egen `cellid_iup' = group( `indid'  `pastunitid' )

quietly tab `timeid'
local periods = r(r)
if (`periods' != 2) {
	display "This specification requires two time periods."
	exit
}

summarize `indid', meanonly 
local nind = r(max) 

summarize `unitid', meanonly 
local nunit = r(max) 

summarize `pastunitid', meanonly 
local npastunit = r(max) 

summarize `cellid_iu', meanonly 
local ncells_iu = r(max) 

summarize `cellid_uup', meanonly 
local ncells_uup = r(max) 

summarize `cellid_iup', meanonly 
local ncells_iup = r(max) 

sort `indid' `unitid'
}
* Launches the Conjugate gradient estimation
* Note : the program fails when explanatory variables are not linearly independent

display " Dependent variable 	: `dependent' "
display " Covariates		: `varlist' "
display " Individual ID 	: `individual' "
display " Unit ID 		: `unit' "
display " Time 			: `time' "

mata: maincgp("`dependent'",tokens("`varlist'"),"`indid'","`unitid'", "`pastunitid'", `nind', `nunit', `npastunit', `ncells_iu', `ncells_uup', `ncells_iup',"`save'")

quietly {

mata: addresultscgp("`dependent'",tokens("`varlist'"),"`indid'","`unitid'", "`pastunitid'", `nind', `nunit', `npastunit', "`save'","`prefix'")

gen resid = y - xb - indeffect - uniteffect - pastuniteffect if pastuniteffect != .
replace resid = y - xb - indeffect - uniteffect  if pastuniteffect == .

summarize indeffect, meanonly
scalar mean_ie = r(mean)
scalar stdev_ie = r(sd)

summarize uniteffect, meanonly
scalar mean_ue = r(mean)
scalar stdev_ue = r(sd)

summarize pastuniteffect, meanonly
scalar mean_upe = r(mean)
scalar stdev_upe = r(sd)

summarize xb, meanonly
scalar mean_xb = r(mean)
scalar stdev_xb = r(sd)

summarize resid, meanonly
scalar mean_resid = r(mean)
scalar stdev_resid = r(sd)

replace indeffect = indeffect - mean_ie

matrix means = (mean_ie, 	mean_ue, 	mean_upe, 	mean_xb, 	mean_resid )
matrix colnames means = Individual_Effect Unit_Effect Past_Unit_Effect XB Residual

matrix sdevs = (stdev_ie, 	stdev_ue, 	stdev_upe, 	stdev_xb, 	stdev_resid)
matrix colnames sdevs = Individual_Effect Unit_Effect Past_Unit_Effect XB Residual

matrix accum R = y indeffect uniteffect resid, nocons dev
matrix R = corr(R)
matrix colnames R = Dependent Individual_Effect Unit_Effect Residual
matrix rownames R = Dependent Individual_Effect Unit_Effect Residual

ereturn scalar nobs 	= _N
ereturn scalar nind     = `nind'
ereturn scalar nunit	= `nunit'
ereturn scalar constant	= mean_ie

ereturn matrix beta	betas
ereturn matrix M		means
ereturn matrix SD		sdevs
ereturn matrix corrtable R
}
end

mata

function maincgp(string scalar dependent, string rowvector covariates,
			string scalar individualid, string scalar unitid, string scalar pastunitid,
			real scalar nind, real scalar nunit, real scalar npastunit,
			real scalar ncells_iu, real scalar ncells_uum, real scalar ncells_ium,
			string scalar filerawoutput) {

///////////////////////////////////////////////////////////////////////
// 
//	STEPS OF THE CONJUGATE GRADIENT ESTIMATION
//    ==========================================
//    
//    1. Build matrices or representations of the matrices :
//		X'X X'D X'F X'F-1 (brutal multiplication by cross)
//		D'D F'F F-1'F-1' (three loops)
//		D'F D'F-1 F'F-1  (building cells - sparse representation)
//		and X'Y (xtprod)
//
//	2. Precondition:
//		Perform the Cholesky decomposition of X'X = U'U
//		Transform matrices off the diagonal
//
//	3. Start conjugate gradient iterations
//		This requires knowing how to multiply by A (modcg)
//
////////////////////////////////////////////////////////////////////////

//*************** DATA DECLARATION

	real scalar n, ncov, ncoef;
	
	real scalar i, j, jp, jf, jfm; 		// Used for storing person and firm indices

	real vector y;				// The vector of observations (Y)					-	Size : n
	real matrix cov;				// The covariate matrix X (view)					- 	Size : n*ncov
	real vector individual;			// The vector of individual indices i				- 	Size : n
	real vector unit;				// The vector of unit indices J(i,t)				- 	Size : n
	real vector pastunit;			// The vector of past unit indices J(i,t-1)			-	Size : n

	real vector d,   f,    fm;          // Sparse representation of matrix D'D F'F F'F-1		-	Size : nind + nunit-1 + npastunit-1

	real vector df,  dfp , dff;        	// Sparse representation of matrix D'F				-	Size : 3*ncell_iu
	real vector dfm, dfmp, dfmfm;      	// Sparse representation of matrix D'F-1				- 	Size : 3*ncell_uum
	real vector ffm, ffmf, ffmfm;      	// Sparse representation of matrix F'F-1				- 	Size : 3*ncell_ium
		
	real matrix xx, xd, xf, xfm;	     	// Matrices X'X X'D X'F X'F-1						-	Size : ncov²+ncov*nind+ncov*(nunit-1)+ncov*(npastunit-1)

	real matrix u ;				// Cholesky decomposition of X'X = U'U				- 	Size : ncov²/2
	real matrix upinv;			// inverse of U'

	real vector b ; 				// Right hand side of the normal form equation of the OLS	-	Size : n

	real vector theta ; 			// The parameter vector 						- 	Size : ncoef = ncov + nind + nunit-1 + npastunit-1

	real scalar icell_iu , icell_uum, icell_ium; 	// Index used to scan through the cells of the dataset

	real scalar maxit; 			// Maximum number of conjugate gradient iterations
	real scalar resid;			// Residual

	real scalar fp;				// File descriptor to save the parameter vector

	real scalar fmissing;			// Value of missing firm	
	
	

//***************** GETTING VIEWS ON THE DATA

	st_view(cov=., 		., covariates);
	st_view(y=., 		., dependent);
	st_view(individual=., 	., individualid);
	st_view(unit=., 		., unitid);
	st_view(pastunit=., 	., pastunitid);

	fmissing = missingof(pastunit);

	n	= length(y);
	ncov  = cols(cov);
	ncoef = ncov + nind + nunit -1 + npastunit -1 ;

	printf("%f observations, %f covariates, %f individuals, %f units, %f past units\n",n,ncov,nind,nunit,npastunit);
	printf("A total of %f coefficients\n",ncoef);	

	theta = J( ncoef, 1, 0);

	xd	= J( ncov,  	nind, 		0);
	xf	= J( ncov, 		nunit-1, 		0);
	xfm	= J( ncov,		npastunit-1, 	0);

	d	= J( nind,        1,    0);
	f 	= J( nunit-1, 	1,    0);
	fm    = J( npastunit-1, 1,    0);

	dfp	= J( ncells_iu,   1,    0);
	dff   = J( ncells_iu,   1,    0);
	df    = J( ncells_iu,   1,    0);

	ffmf	= J( ncells_uum,   1,    0);
	ffmfm = J( ncells_uum,   1,    0);
	ffm   = J( ncells_uum,   1,    0);

	dfmp	= J( ncells_ium,   1,    0);
	dfmfm = J( ncells_ium,   1,    0);
	dfm   = J( ncells_ium,   1,    0);

	b	= J( ncoef,		 1,	 0);



	if ((length(y) != rows(cov)) | (length(y)!=length(individual)) | (length(y)!=length(unit)) | (length(y)!=length(pastunit)) ) {
		printf("Error : invalid matrix length\n");
		exit();
	}


	xx = cross(cov,cov);

	for ( i=1; i<=n ; i++ ) {
		jp 		= individual[i]
		d[jp] 	= d[jp] + 1
		
		jf 		= unit[i]
		if (jf != nunit) {
			f[jf] 	= f[jf] + 1 
			assert(jf < nunit);
		}

		jfm 		= pastunit[i]
		if (jfm != fmissing && jfm >= 1 && jfm != npastunit) {
			fm[jfm] 	= fm[jfm] + 1
			assert(pastunit[i] < npastunit);
		}

		for ( j = 1; j <= ncov ; j++ ) {
			/* Construction de XD XF et XF-1 					*/
			/* Cela suppose que xd xf et xfm soient initialisés à 0 	*/
			xd[j,jp]   = xd[j,jp] 	+ cov[i, j]
			if (jf != nunit) {
				xf[j,jf]   = xf[j,jf] 	+ cov[i, j];
			}
			if (jfm != fmissing && jfm >= 1 && jfm != npastunit) {
				xfm[j,jfm] = xfm[j,jfm] + cov[i, j];
			}
		}

	}

//************ BUILDING SPARSE REPRESENTATIONS OF D'F F'F-1 D'F-1

// D'F

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
		}
	}

	assert(icell_uum == ncells_uum);


//********** PRECONDITIONING

	printf("Starting preconditioning\n");

	u = cholesky(xx)';
	upinv = luinv(u');

	cov	= cov * upinv ;

	xd 	= upinv * xd;
	xf 	= upinv * xf;
	xfm 	= upinv * xfm ;

	for (i = 1; i<=nind ; i ++) {
		d[i] = 1/sqrt(d[i]);
	}

	for (i = 1; i<=nunit-1 ; i ++) {
		f[i] = 1/sqrt(f[i]);
	}

	for (i = 1; i<=npastunit-1 ; i ++) {
		fm[i] = 1/sqrt(fm[i]);
	}


	for ( i = 1 ; i <= ncov ; i++) {
		for (j = 1; j <= nind ; j ++ ) {
			xd[i, j] = xd[i, j]	*	d[j];
		}
		for (j = 1; j <= nunit-1 ; j ++ ) {
			xf[i, j] = xf[i, j]	*	f[j];
		}
		for (j = 1; j <= npastunit-1 ; j ++ ) {
			xfm[i, j] = xfm[i, j]	*	fm[j];
		}

	}


	for (i=1; i<=ncells_iu; i++) {
		if (dff[i] < nunit) {	
			df[i] = df[i]*d[dfp[i]]*f[dff[i]]
		}
	}

	for (i=1; i<=ncells_ium; i++) {
		if (dfmfm[i]<npastunit) {
			dfm[i] = dfm[i]*d[dfmp[i]]*fm[dfmfm[i]]
		}
	}

	for (i=1; i<=ncells_uum; i++) {	
		if (ffmf[i] < nunit && ffmfm[i] < npastunit ) {
			ffm[i] = ffm[i]*f[ffmf[i]]*fm[ffmfm[i]]
		}
	}

	xtprod (individual, unit, pastunit, cov, y, d, f, fm, b, 
			n, ncov, nind, nunit, npastunit );

	resid = 1.0e-7
	maxit = 1000


	printf("Beginning Conjugate Gradient Iterations\n");


	modcg(b, theta, maxit, resid, xd, xf, xfm, 
		df, dfp, dff, ncells_iu,
		dfm, dfmp, dfmfm,	ncells_ium,
		ffm, ffmf, ffmfm,	ncells_uum,
		ncov, nind, nunit, npastunit, ncoef, n);

	printf("Transforming theta back to original scale\n");


	theta[1 .. ncov] = luinv(u) * theta[1 .. ncov];

	for (jp = 1; jp <= nind; jp++) {
		j = ncov + jp
		theta[j] = theta[j] * d[jp]
	}


	for (jf = 1; jf < nunit; jf++) {
		j = ncov + nind + jf
		theta[j] = theta[j] * f[jf]
	}

	for (jfm = 1; jfm < npastunit ; jfm++) {
		j = ncov + nind + nunit-1 + jfm
		theta[j] = theta[j] * fm[jfm]
	}

	//printf("Saving to file %s\n",filerawoutput);
	fp = fopen(filerawoutput,"w");

	fputmatrix(fp, theta);

	fclose(fp);


}


void xtprod(real vector individual, real vector unit, real vector pastunit, real matrix cov, // the X = (beta theta psi psi-1) matrix
		real vector y, real vector d, real vector f, real vector fm, real vector r,
		real scalar n, real scalar ncov, real scalar nind, real scalar nunit, real scalar npastunit
		)

///////////////////////////////////////////////////////////
//
// MULTIPLIES X's = r
//
///////////////////////////////////////////////////////////

{

//////////////////////////
//		( beta's  ) //
// X's  = 	( theta's ) //
//		( psi's   ) //
//		( psi-1's ) //
//////////////////////////

	real scalar i, jind, junit, jpastunit, tmp, j, k, l;	
	real scalar fmissing;
	real scalar ncoef;

	fmissing = missingof(pastunit);

	ncoef = ncov + nind + nunit - 1 + npastunit - 1 ;
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
		l = ncov + nind + nunit - 1 + jpastunit 

		
		r[ j ] = r[ j ] +	d[jind] 	  * tmp;

		if ( junit <= nunit - 1) {
			r[ k ] = r[ k ] +	f[junit] 	  * tmp;
		}

		if ( jpastunit != fmissing && jpastunit <=npastunit - 1) {   // Only if the past unit id is less than the maximum and not missing
			r[ l ]  = r[ l ]    + fm[jpastunit] * tmp;
		}
	}

}

void modcg(	real vector b, real vector x, real scalar maxit, real scalar resid, 
		real matrix xd, real matrix xf, 	real matrix xfm,
		real matrix df, real matrix dfp, 	real matrix dff,	real scalar ncells_iu,
		real matrix dfm, real matrix dfmp, real matrix dfmfm,	real scalar ncells_ium,
		real matrix ffm, real matrix ffmf, real matrix ffmfm,	real scalar ncells_uum,
		real scalar ncov, real scalar nind, real scalar nunit, real scalar npastunit,
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
	
	matvec (	x, r, 
			xd, 	xf,	xfm,
			df, 	dfp,	dff, 		ncells_iu,
			dfm,	dfmp,	dfmfm,	ncells_ium,
			ffm,  ffmf,	ffmfm,	ncells_uum,
			ncov, nind, nunit, npastunit,	ncoef );

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
		matvec (p, q, xd, 	xf,	xfm,
			df, 	dfp,	dff, 		ncells_iu,
			dfm,	dfmp,	dfmfm,	ncells_ium,
			ffm,  ffmf,	ffmfm,	ncells_uum,
			ncov, nind, nunit, npastunit,	ncoef);
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


void matvec (
		real vector xin, real vector rout,
		real matrix xd,   real matrix xf,	real matrix xfm,					
		real matrix df,   real matrix dfp,	real matrix dff,		real scalar ncells_iu,
		real matrix dfm,  real matrix dfmp,	real matrix dfmfm,	real scalar ncells_ium,
		real matrix ffm,  real matrix ffmf,	real matrix ffmfm,	real scalar ncells_uum,
		real scalar ncov, 
		real scalar nind, 
		real scalar nunit,
		real scalar npastunit,
		real scalar ncoef
		)

/////////////////////////////////////////////////////
//                                                 //
// MATVEC computes rout = A*xin                    //
//                                                 //
// with A = ( X D F F-1 )' ( X D F F-1 )           //
//								   //
/////////////////////////////////////////////////////

{
	
	real scalar i, jind, junit, jpastunit;

///////////////////////////////////////////////////////////////////
//		( beta  + X'D theta  + X'F psi     + X'F-1 psi-1 )   //
//  AX =    ( theta + D'X beta   + D'F psi     + D'F-1 psi-1 )   //
//		( psi   + F'X beta   + F'D theta   + F'F-1 psi-1 )   //
//		( psi-1 + F-1'X beta + F-1'D theta + F-1'F psi   )   //
// 										     //
//	X'D, X'F, X'F-1 are stored as real matrices.               //
//    Other matrices are stored as sparse matrices.		     //
///////////////////////////////////////////////////////////////////

	rout = xin;

// First line
	rout[1 .. ncov] = xd*xin[ncov+1 .. ncov+nind] 								+ rout[1 .. ncov]
	rout[1 .. ncov] = xf*xin[ncov+nind+1 .. ncov+nind+nunit-1] 						+ rout[1 .. ncov]
	rout[1 .. ncov] = xfm*xin[ncov+nind+(nunit-1)+1 .. ncov+nind+(nunit-1)+npastunit-1]  	+ rout[1 .. ncov]

// Second line without sparse matrices : adds X'D beta
	rout[ncov+1 .. ncov+nind ] = cross(xd, xin[1 .. ncov]) 			+ rout[ncov+1 .. ncov+nind] 

// Third line without sparse matrices  : adds F'X beta
	rout[ncov+nind+1 .. ncov+nind+nunit-1] = cross(xf, xin[1 .. ncov]) 	+ rout[ncov+nind+1 .. ncov+nind+nunit - 1]

// Fourth line without sparse matrices : adds F-1'X beta
	rout[ncov+nind+(nunit - 1)+1 .. ncov + nind+(nunit-1)+npastunit-1] = cross(xfm, xin[1 .. ncov]) + rout[ncov+nind+(nunit-1)+1 .. ncov+nind+(nunit-1)+npastunit-1] 

// Adding the contents of D'F
	for (i = 1; i<= ncells_iu; i++) {
		jind  = ncov + dfp[i];
		junit = ncov + nind   + dff[i];

		if (dff[i] != nunit) {
			rout[jind]  = xin[junit]*df[i] + rout[jind]
			rout[junit] = xin[jind] *df[i] + rout[junit]
		}
	}

// Adding the contents of D'F-1
	for (i = 1 ; i <= ncells_ium ; i++ ) {
		jind 		= ncov + dfmp[i];
		jpastunit 	= ncov + nind + nunit - 1 + dfmfm[i]
		
		if (dfmfm[i] != npastunit) {
			rout[jind] 		= rout[jind] + xin[jpastunit] * dfm[i]
			rout[jpastunit] 	= rout[jpastunit] + xin[jind] * dfm[i]
		}
	} 
	
// Adding the contents of F'F-1
	for (i = 1 ; i <= ncells_uum ; i++ ) {
		junit		= ncov + nind + ffmf[i] ;
		jpastunit	= ncov + nind + nunit - 1 + ffmfm[i] ;

		if (ffmf[i] != nunit && ffmfm[i] != npastunit) {
			rout[junit]		= rout[junit] + xin[jpastunit] * ffm[i]
			rout[jpastunit]	= rout[jpastunit] + xin[junit] * ffm[i]
		}
	}
}

function addresultscgp(string scalar dependent, string rowvector covariates,
				string scalar individualid, string scalar unitid, string scalar pastunitid,
				real scalar nind, real scalar nunit, real scalar npastunit,
				string scalar filerawoutput,
				string scalar prefix) {


	real scalar ncov;
	real scalar fp;
	real vector betas, indeffect, uniteffect, pastuniteffect;
	real vector params;
	real matrix data;
	real scalar n,i;
	real vector indeffect_name, uniteffect_name, pastuniteffect_name;
	string scalar cmd_listcovnames;
	string scalar cmd_listcovvalues;

	ncov = length(covariates);	

	fp = fopen(filerawoutput,"r")
	params = fgetmatrix(fp);
	
	betas 		= params[1 .. ncov];
	indeffect 		= params[ncov+1 .. ncov + nind];
	uniteffect 		= params[ncov+nind+1 .. ncov + nind + nunit -1];
	pastuniteffect 	= params[ncov+nind+(nunit-1)+1 .. ncov + nind + (nunit -1) + (npastunit - 1)];

	indeffect_name = "indeffect";
	uniteffect_name = "uniteffect";
	pastuniteffect_name = "pastuniteffect";
	st_addvar("double", indeffect_name);
	st_addvar("double", uniteffect_name);	
	st_addvar("double", pastuniteffect_name);
	st_view(data, . ,(individualid, unitid, pastunitid, indeffect_name, uniteffect_name, pastuniteffect_name));

	n = rows(data);

	for (i = 1; i<= n ; i ++) {
		data[i,4] = indeffect[data[i,1]];
		if (data[i,2] != nunit) {
			data[i,5] = uniteffect[data[i,2]];
		} else {
			data[i,5] = 0;
		}
		if (data[i,3] != npastunit && data[i,3] != .) {
			data[i,6] = pastuniteffect[data[i,3]];
		} else if (data[i,3] != .) {
			data[i,6] = 0;
		} else {
			data[i,6] = .;
		}
	}

	stata("gen xb = 0");

	cmd_listcovnames 	=	"";
	cmd_listcovvalues = 	"";

	for (i = 1; i< ncov; i++) {
		cmd_listcovnames 	= sprintf("%s %s",cmd_listcovnames, covariates[i]);
		cmd_listcovvalues	= sprintf("%s %g,",cmd_listcovvalues, betas[i]);
		stata(sprintf("replace xb = xb + %g * %s",betas[i], covariates[i]));
	}

	cmd_listcovnames 	= sprintf("%s %s",cmd_listcovnames, covariates[ncov]);
	cmd_listcovvalues	= sprintf("%s %g",cmd_listcovvalues, betas[ncov]);
	stata(sprintf("replace xb = xb + %g * %s",betas[ncov], covariates[ncov]));

	stata(sprintf("matrix input betas = (%s)",cmd_listcovvalues));
	stata(sprintf("matrix colnames betas = %s",cmd_listcovnames));
	
}


end

