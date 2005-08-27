*** Pull the necessary variables from the annearn06 file ***;
*** Create the binary file ***;

%include './config_param.sas';

options obs=max mprint mlogic symbolgen;

%bincalc(4);

/* Not necessary to create the basefile     */
/* data cg.basefile() cg.year(keep=year);   */
/* To create the basefile use the line above*/

*data _NULL_ cg.year(keep=year);
 data cg.basefile() cg.year(keep=year); 
 
  merge cg.wage_history_01 cg.strip02;
     by pik sein year;

  /* Can place other data manipulation statements here */
     
  /* Output the binary file on the fly */

  file "&cellout./cgout4" lrecl=&binlen recfm=N;

     put (&binvars.) (float4.)
         (persseq firmseq) (ib4.);
     output cg.basefile;
     output cg.year;
run;

/* No means check since this is handled by CG2       */
/* and for large problems this is costly             */
/* Uncomment to check if CG2 read in data correctly  */
   
*proc contents data=cg.basefile;
*proc means data=cg.basefile;
*  output out=dot.means04 mean=;
*run;
   
/* Data Quality Checks */

*proc print data=cg.basefile(obs=100);
*  id pik sein year;
*run;
