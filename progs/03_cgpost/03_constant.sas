*** Calculate the constant ***;

%include './config_param.sas';


data mean0;
   set dot.means02(keep=&rhs);
run;

proc print;
run;

data beta0;
   set dot.rhs(keep=&rhs);
run;

proc print;
run;

data y0;
   set dot.means02(keep=&depvar);
run;

proc print;
run;


proc iml;

 /*Read in xbar*/
 use mean0 var {&rhs};
   read current into xbar;
 close;
 /*Read in the betas*/
 use beta0 var {&rhs};
   read current into beta;
 close;
 /*Read in ybar */
 use y0 var _all_;
   read current into ybar;
 close;

 print xbar beta ybar;

 xb=xbar*beta`;

 print xb;

 cons=ybar-xb;

 print cons;

 create dot.constant var{cons};
   append;
   show contents;
 close dot.constant;
quit;
run;
