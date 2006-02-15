*** Create necessary environment for post processing CG estimates **;

%include './config_param.sas';

/* Create the cg.coef file */

x "wc -c config_param.sas | awk ' { print $1 }' > 00_setup.dat";

data _NULL_;
  infile '00_setup.dat';

  input maxlen;

  call symput('maxlen',maxlen);
run;

x 'rm -f 00_setup.dat';

%put maxlen=&maxlen;

data _NULL_;

file 'cg.coef' lrecl=&maxlen;

put "&depvar" " " "&persid" " " "&firmid" " " "&rhs";
run;

/* Create a link to the betas file */

x 'rm -f cg.betas';
x 'rm -f cg.cgin';
x "ln -s &betadir./cg.betas";
x "ln -s &betadir./cg.cgin";

/* Create a link to the means file */
/* This file only exists for the v3 version */
/* The link will be created in all cases even */
/* if the file does not exist */

x 'rm -f cg.means';
x "ln -s &betadir./cg.means";
