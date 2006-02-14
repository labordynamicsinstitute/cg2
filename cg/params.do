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
global personid ="pupilid"
global firmid   ="schoolid"
global CGIN ="cgin"
global workingdir = "W:\projets\schooleffects\src\cg"

global DEPENDENT = "y"
global COVARIATES = "keystage1 keystage2"
*** This should really be computed 
global NCOV = 1

cd $workingdir

********* DO NOT MODIFY AFTER THIS LINE

use $CGIN, clear
keep $personid
duplicates drop $personid, force
global NPUPILS = "`=_N'"

use $CGIN, clear
keep $firmid
duplicates drop $firmid, force
global NSCHOOLS = "`=_N'"

use $CGIN, clear
keep $personid $firmid
duplicates drop $personid $firmid, force
global NCELLS = "`=_N'"

clear
program drop _all

do "$VARFILE.do"

defvar
