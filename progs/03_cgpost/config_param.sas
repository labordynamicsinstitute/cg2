
options obs=max symbolgen mprint mlogic;

libname dot '.';

*** Macro Variables Necessary to create the HC File ***;

%let depvar=learn;
%let persid=pik;
%let firmid=sein;

%let betadir=/home/kmckinney/runcg8_out;

%let cgin=cgout4;

/* Only 4 byte and 8 byte are acceptable*/
%let cginlen=4;

%let BSTART=1;

libname cgdir "&betadir";

/* This is the same as binvars except depvar is not included */

%let rhs=exper;

*** Create several important macro variables using cg.cgin ***;
** NN     = Number of Observations **;
** NCELLS = Number of Cells of Person Firm Matches **;
** NNPERS = Number of Persons **;
** NNFIRM = Number of Firms **;
** NCOV1  = Number of RHS vars **;
************************************************;

data _NULL_;
  infile "./cg.cgin";

  input nn ncells nnpers nnfirm ncov1;

  call symput('NN',nn);
  call symput('NCELLS',ncells);
  call symput('NNPERS',nnpers);
  call symput('NNFIRM',nnfirm);
  call symput('NCOV1',ncov1);
run;

%put _ALL_;

** This version is slightly different than the versions in runcg ***;
*** Here I calculate the length of a record in the binary file ***;
*** However the number of vars is only the RHS vars with the depvar assumed ***;

%macro bincalc(len);
   %global binlen;
   %let i=1;
   %let v=%scan(&rhs,&i);
   %do %until ("&v"="");
     %let i=%eval(&i+1);
     %let v=%scan(&rhs,&i);
   %end;
   %let binlen=%eval(8+&len.*(&i.-1)+&len);
   %put "bin length= " &binlen;
   %put "num coef= " %eval(&i-1);
%mend;

%macro rhs_tmp; 
   %global rhs_tmp num_rhs;
   %let i=1;
   %let v=%scan(&rhs,&i);
   %let rhs_tmp=&v._tmp;
   %do %until ("&v"="");
     %let i=%eval(&i+1);
     %let v=%scan(&rhs,&i);
     %if ("&v"~="") %then %do;
        %let rhs_tmp=&rhs_tmp &v._tmp;
     %end;
   %end;
   %put "rhs_tmp=" &rhs_tmp;
   %let num_rhs=%eval(&i-1);
%mend;

/* read in betas */
%macro vars ;
   %let i=&BSTART;
   %do %until ("&variable"="");
       %let variable=%scan(&rhs,&i,' ');
       %global &variable  ;
       %if ("&variable"~="") %then %do;
          call symput("&variable",&variable);
       %end;
       %let i=%eval(&i+1);
   %end;
%mend;

/* calculate xb without the constant */
%macro calc_xb;
%let i=&BSTART;
   %do %until ("&variable"="" ) ;
       %let variable=%scan(&rhs,&i,' ');
       %if ("&variable"~="") %then %do;
          + &variable * &&&variable
       %end;
       %let i=%eval(&i+1);
   %end;
%mend;

/* Create the correct input statement for 4 and 8 byte data */
%macro inputbin;
      %if &cginlen=4 %then %do;
         input (&depvar._tmp &rhs_tmp) (float4.)
            (persseq firmseq) (ib4.);
      %end;
      %else %do;
         input (&depvar._tmp &rhs_tmp) (rb8.)
            (persseq firmseq) (ib4.);
      %end;
%mend;

/* Create the correct input statement for 4 and 8 byte data */
%macro inbin;
      %if &cginlen=4 %then %do;
         input (&depvar &rhs) (float4.)
            (persseq firmseq) (ib4.);
      %end;
      %else %do;
         input (&depvar &rhs) (rb8.)
            (persseq firmseq) (ib4.);
      %end;
%mend;

