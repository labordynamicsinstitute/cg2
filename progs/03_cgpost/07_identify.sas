**** Identify the model and fix grouping problems ****;

%include './config_param.sas';

proc sort data=cgdir.hcest1 out=temp1;
   by group;
run;

proc standard data=temp1 out=temp1 mean=0;
   by group;
   var theta;
run;

proc standard data=temp1 out=temp1 mean=0;
   var psi;
run;

proc sort data=temp1;
   by pik sein;
run;

proc means data=temp1(where=(group_size=1)) n mean stddev ;
   title1 "Stats on theta psi before fix";
   var theta psi;
run;

/* Fix the unidentified person and firm effects */
/* When groupsize=1 randomly assign a person    */
/* and firm effect using a normal approximation */
/* to the actual theta and psi distribution     */

proc means data=temp1 n mean stddev;
   title1 "Stats on theta psi for all workers and firms";
   var theta psi;
   output out=sumstats mean(theta psi)=mean_theta mean_psi
                       stddev(theta psi)=sd_theta sd_psi;
run;

proc print data=sumstats;
run;
   
/* Output the final file */

data cgdir.hcest2(drop=mean_theta mean_psi sd_theta sd_psi draw_theta draw_psi);
   set temp1;
      by pik;

   retain draw_theta draw_psi;
   if _N_=1 then set sumstats(keep=mean_theta mean_psi sd_theta sd_psi);

   if first.pik then do;
      draw_theta=ranuni(15389);
      draw_psi=ranuni(853762);
   end;
   if group_size=1 then do;
      theta=mean_theta+probit(draw_theta)*sd_theta;
      psi=mean_psi+probit(draw_psi)*sd_psi;
   end;

   h=cons+exper+theta;
   resid=&depvar-cons-theta-psi-xb;
run;

proc means data=cgdir.hcest2(where=(group_size=1)) n mean stddev ;
   title1 "stats on theta and psi after fix";
   var theta psi;
run;

proc contents;
   title1 "Summary Stats and Correlations for all workers";
proc corr data=cgdir.hcest2;
   var &depvar theta psi xb exper h resid;
proc print data=cgdir.hcest2(obs=100);
  id pik sein;
run;
