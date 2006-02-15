*** Config file for creation of CG input file ***;

** binvars=depvar rhs **;
** perseq and firmseq are automatically included **;

%let binvars=learn exper;

%put "binvars=" &binvars;

** Location of input data files **;

   libname inputs '/home/kmckinney/ciser_data';

*** location of the output directory ***;

   %let cellout=/home/kmckinney/runcg8_out;

*** Output ***;

libname cg "&cellout.";
libname dot '.';

*** MACROS ***;

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

*** Calculate the length of a record in the binary file ***;

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
