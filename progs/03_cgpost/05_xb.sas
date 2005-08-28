*** Create the xb index for each record without the constant ***;

%include './config_param.sas';

%bincalc(&cginlen.);

/*Put each beta value into a macro variable*/

data _null_;
  set dot.rhs;
  %vars ; 
run;

%put _GLOBAL_;

data xb(keep=persseq firmseq &depvar xb exper);

   infile "&betadir./&cgin." lrecl=&binlen recfm=N;

   %inbin;

   xb=0 %calc_xb;

   exper=exper * &exper;

 /* Usually we use a quartic in experience */
 /*        exper2* &exper2+
         exper3* &exper3+
         exper4* &exper4; */
run;

/* Attach the year to each record */

data cgdir.xb;
   merge xb cgdir.year;
run;

proc print data=cgdir.xb(obs=1000);
run;
