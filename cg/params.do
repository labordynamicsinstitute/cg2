
***** Sets parameters for estimation
*****
*****

** Change line here if you want to include different covariates
** !! Be careful no check is performed on the consistency between
** !! and the string

global VARFILE ="variables"

cd "W:\projets\schooleffects\src\cg"

********* DO NOT MODIFY AFTER THIS LINE

use cgin, clear
keep pupilid
duplicates drop pupilid, force
global NPUPILS = "`=_N'"

use cgin, clear
keep schoolid
duplicates drop schoolid, force
global NSCHOOLS = "`=_N'"

use cgin, clear
keep pupilid schoolid
duplicates drop pupilid schoolid, force
global NCELLS = "`=_N'"

clear
program drop _all

do "$VARFILE.do"

defvar
