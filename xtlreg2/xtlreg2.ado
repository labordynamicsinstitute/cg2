cap program drop xtlreg2
program define xtlreg2, eclass
* XTLREG uses mata, therefore requires version 9
version 9

global 	AUTHOR	=	"Amine Ouazad"
global	VERSION	=	"0.1"

/* 
 * 
 * Syntax: xtlreg2 depvar indepvars, individual(first fixed effect) 
 *		unit(second fixed effect) pastunit(third fixed effect) lambda(value of the dynamic coefficient) 
 *
 */

syntax varlist(min=2), individual(varname) unit(varname) lambda(real) time(varname) sorted [prefix(name)] [save(name)]

display "xtlreg2 program version $VERSION - Estimation of two-way fixed effects with past unit effects & past individual effects"

cap drop indeffect uniteffect pastuniteffect xb resid pred

/* Checks whether all variables are numeric */
display "Checks whether all variables are numeric"

foreach v in `varlist' {
	confirm numeric variable `v'
}

gettoken dependent varlist : varlist

tempvar indid unitid pastunitid timeid cellid_iu cellid_uup cellid_iup t

if ("`save'"=="") {
	tempfile save
}

if ("`prefix'"=="") {
	local prefix "effect"
}

* Creating sequenced ID variables
display "Creating sequenced ID variables"

egen `indid'  = group(`individual')
egen `unitid' = group(`unit')
egen `timeid' = group(`time')

quietly {
	tab `timeid'
	local number_years = r(r)
	if (`number_years' != 2) {
		error "There should be two periods in the dataset"
	}
}

tempfile temp_ids temp_data
save "`temp_data'", replace

keep `indid' `timeid' `unitid'
duplicates drop `indid' `timeid', force
tsset `indid' `timeid'
gen `pastunitid' =l.`unitid'
sort `indid' `timeid'
save "`temp_ids'", replace

use "`temp_data'"
sort `indid' `timeid'
merge `indid' `timeid' using "`temp_ids'", keep(`indid' `timeid' `pastunitid')
drop _merge

egen `cellid_iu' = group(`indid' `unitid')

sort `unitid' `pastunitid'
egen `cellid_uup' = group( `unitid' `pastunitid' ) 

sort `indid' `pastunitid'
egen `cellid_iup' = group( `indid'  `pastunitid' )

summarize `indid', meanonly 
local nind = r(max) 

summarize `unitid', meanonly 
local nunit = r(max) 

summarize `cellid_iu', meanonly 
local ncells_iu = r(max) 

summarize `cellid_uup', meanonly 
local ncells_uup = r(max) 

summarize `cellid_iup', meanonly 
local ncells_iup = r(max) 

sort `indid' `unitid'
*}
* Launches the Conjugate gradient estimation
* Note : the program fails when explanatory variables are not linearly independent

display " Dependent variable 	: `dependent' "
display " Covariates		: `varlist' "
display " Individual ID 	: `individual' "
display " Unit ID 		: `unit' "
display " Time 			: `time' "

mata: maincgxtl2("`dependent'",tokens("`varlist'"),"`indid'","`unitid'", "`pastunitid'", "`timeid'", `lambda', `nind', `nunit', `ncells_iu', `ncells_uup', `ncells_iup',"`save'")

quietly {

mata: addresultscgxtl2("`dependent'",tokens("`varlist'"),"`indid'","`unitid'", "`pastunitid'", `nind', `nunit', "`save'","`prefix'")

gen pred = xb + (1+`lambda')*indeffect +  uniteffect + `lambda' * pastuniteffect if `timeid' == 2
replace pred = xb + indeffect + uniteffect  if `timeid' == 1

gen resid = y - pred

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

scalar RSS = mean_resid * _N

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

