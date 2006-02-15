*** $Id$ ;
*** $URL$ ;
***;
*** Original author: Kevin McKinney ;
***;
*** Pull the necessary variables from the annearn06 file ***;
*** Create the binary file ***;

%include './config_param.sas';

options obs=max mprint mlogic symbolgen;

/*------------------------------------------------------------
    All parameters can be adjusted at the end of this file
------------------------------------------------------------*/

/*BEGINCCC
   Calculate the length of a record in the binary file for a 
   given width of each variable. Result is stored in a global 
   macro variable 'binlen'.
CCCEND*/

%bincalc(4);

/*BEGINCCC
	Local macro to run things.
CCCEND*/
/*------------------------------------------------------------
 Syntax: (defaults in CAPS)

  %create_binary(
              diagnostics=yes/NO,   Detailed diagnostics
              basefile=yes,NO       Permanent basefile (not needed)
	      )
------------------------------------------------------------*/

%macro create_binary(diagnostics=no,basefile=no);

/*BEGINCCC
  It is not necessary to create the basefile  
  To actually create the basefile use the line below:
CCCEND*/
  
%if ( &basefile = yes ) %then %do;
       data MYWORK.basefile() MYWORK.year(keep=year);   
%end;
%else %do;
       data _NULL_ MYWORK.year(keep=year);
%end;

  merge MYWORK.wage_history_01 MYWORK.strip02;
     by pik sein year;

/*BEGINCCC
   Can place other data manipulation statements here.
CCCEND*/

/*   The macro is defined in site_macros.sas */

     %site_manipulations;

/*    Output the binary file on the fly */

  file "&cellout./cgout4" lrecl=&binlen recfm=N;

     put (&binvars.) (float4.)
         (persseq firmseq) (ib4.);
     output MYWORK.basefile;
     output MYWORK.year;
run;

/* No means check since this is handled by CG2       */
/* and for large problems this is costly             */

%if ( &diagnostics = yes ) %then %do;   
   %put RUNNING DIAGNOSTICS;
   proc contents data=MYWORK.basefile;
   proc means data=MYWORK.basefile;
   output out=dot.means04 mean=;
   run;
   
/* Data Quality Checks */

   proc print data=MYWORK.basefile(obs=100);
     id pik sein year;
   run;

%end; /* end of diagnostics */
%mend;

/*------------------------------------------------------------
      Now run it 
------------------------------------------------------------*/
%create_binary(diagnostics=no,basefile=no);