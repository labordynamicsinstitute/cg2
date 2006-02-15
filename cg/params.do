* $Id$ 
* $URL$
*****
***** Sets parameters for estimation
*****
***** Requires the file $VARFILE.do
*****

** Change line here if you want to include different covariates
** !! Be careful no check is performed on the consistency between
** !! and the string

global VARFILE ="variables"
global personid ="pik"
global firmid   ="sein"
global CGIN ="wage_history_01"
global workingdir = "."

global DEPENDENT = "learn"
*** This is replicated in variables.do
global COVARIATES = "exper year"
*** This should really be computed 
global NCOV = 2

cd $workingdir

********* DO NOT MODIFY AFTER THIS LINE

use $personid using $CGIN, clear
duplicates drop $personid, force
global NPERSONS = "`=_N'"

use $firmid using $CGIN, clear
duplicates drop $firmid, force
global NFIRMS = "`=_N'"

use $personid $firmid using $CGIN, clear
duplicates drop $personid $firmid, force
global NCELLS = "`=_N'"

clear
program drop _all

do "$VARFILE.do"

di "NPERSONS = $NPERSONS"
di "NFIRMS   = $NFIRMS"
di "NCELLS   = $NCELLS"

defvar
