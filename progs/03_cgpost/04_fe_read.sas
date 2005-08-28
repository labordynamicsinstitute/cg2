*** Read in the group psi and theta for each cell ***;

%include './config_param.sas';


/* read in groups output */
data groups;
   infile "&betadir./groups/groups";
   length persseq firmseq group 5;
   input persseq firmseq group;
run;

proc contents;
proc print data=groups(obs=10);
run;

/* Get the Group Size */
proc sort data=groups;
   by group;
run;

data grptmp(keep=group group_size);
   set groups;
      by group;

   retain group_size;

   if first.group then group_size=0;

   group_size=group_size+1;

   if last.group then output;
run;

data groups;
   merge groups grptmp;
      by group;
run;

proc sort data=groups out=cgdir.groups;
   by persseq firmseq;
run;

/* read in thetas */
%let offset=&NCOV1 ; 
data cgdir.theta(keep=persseq theta) cgdir.psi(keep=firmseq psi) ;
     length firmseq persseq 5;
     infile "&betadir./cg.betas"
           firstobs=%eval(&offset + 1)
         ;
     input var1 var2 ;
     if ( _n_ > %eval( &NNPERS ) ) then do;
         firmseq = var1;
         psi = var2;
         output cgdir.psi;
         end;
     else do;
         persseq = var1;
         theta = var2;
         output cgdir.theta;
         end;
run;

proc print data=cgdir.theta(obs=10);
run;
proc print data=cgdir.psi(obs=10);
run;

/* Run some checks on the data */

proc means data=cgdir.groups N min max;
proc means data=cgdir.theta N min max;
proc means data=cgdir.psi N min max;
run;
