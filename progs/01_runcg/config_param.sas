*** $Id$ ;
*** $URL$ ;
***;
*** Original author: Kevin McKinney ;
***;
*** Config file for creation of CG input file ***;

/*BEGINCCC
	The binvars macro variable should contain a list beginning with the
	dependent variable followed by each right hand side variable (do not
	include a constant, do not include the identifiers).
CCCEND*/

** Syntax: binvars=depvar rhs **;
** perseq and firmseq are automatically included **;

%let binvars=learn exper;

** Location of input data files **;

   libname inputs '/path/to/my/inputs';

*** location of the output from SAS, input to CG ***;

   %let cellout=../02_runcg_out;

*** Working output ***;

   libname mywork ".";
   libname dot ".";

*** MACROS ***;

/*BEGINCCC
	 Calculate the length of a record in the binary file
CCCEND*/

%macro bincalc(len);
   %global binlen;
   %let i=1;
   %let v=%scan(&binvars,&i);
   %do %until ("&v"="");
     %let i=%eval(&i+1);
     %let v=%scan(&binvars,&i);
   %end;
   %let binlen=%eval(8+&len.*(&i.-1));
   %put "bin length= " &binlen;
   %put "num coef= " %eval(&i-1);
%mend;

/*BEGINCCC
	Put all macros specific to your site in 
        this external file.
CCCEND*/

%include "./site_macros.sas";

/*BEGINCCC
	Generate some output.
CCCEND*/

%put "binvars=" &binvars;
