*** $Id$ ;
*** $URL$ ;
***;
*** Original author: Kevin McKinney ;
***;

%macro site_manipulations;

 %put Placeholder for local data manipulations.;
 %put Adjust site_macros.sas to your needs.;

%mend;


*** Create libnames ***;
%macro aelibs;
   %let i=1;
   %let st=%scan(&states,&i);
   %do %until ("&st"="");
     libname ae&st "/data/master/hc_estimates/v1.2/&st";
     libname ehf&st "/data/production/prod/current/ehf/&st";
     libname icf&st "/data/production/prod/current/icf/&st";
     %let i=%eval(&i+1);
     %let st=%scan(&states,&i);
   %end;
%mend;

*** Interleave set statement ***;
%macro set_ae;
   set
   %let i=1;
   %let st=%scan(&states,&i);
   %do %until ("&st"="");
      /* Create the set statement */
      ae&st..annearn06 (keep=&keepae6)
     %let i=%eval(&i+1);
     %let st=%scan(&states,&i);
   %end;;
%mend;

*** Data view ***;

*** Macro to Create a Data View of the Annearn Files ***;

%macro aeview;
   data aestack /view=aestack;
   %let i=1;
   %let st=%scan(&states,&i);
   set
   %do %until ("&st"="");
     ae&st..annearn06(keep=&keepae6)
     %let i=%eval(&i+1);
     %let st=%scan(&states,&i);
   %end;;
     by pik sein year;
   run;
%mend;

