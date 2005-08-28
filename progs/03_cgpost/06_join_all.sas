*** Create the xb index for each record without the constant ***;

%include './config_param.sas';


/* Sort the group file by firmseq */
proc sort data=cgdir.groups(keep=group firmseq persseq group_size)
          out=grptmp;
by firmseq;
run;

/* Atach psi to the group file */
data psi ;
    merge cgdir.psi grptmp;
    by firmseq; 
    if psi=. then psi=0;
run;
proc sort data=psi;
   by persseq firmseq;
run;

/* Attach the PIK and SEIN to the theta file */
data fixed;
     merge cgdir._lookup 
           cgdir.theta 
         ;
     by persseq;
run;

/* Bring the group psi and theta information together */
data fixed;
     merge fixed
           psi
         ;
     by  persseq firmseq;
run;

/* Attach everything to the xb file */;

data cgdir.hcest1;
     merge cgdir.xb fixed;
       by persseq firmseq;

     if _N_=1 then set dot.constant;
run;

proc contents;
proc means data=cgdir.hcest1 N mean min max;
proc print data=cgdir.hcest1(obs=100);
  id pik sein;
print;
