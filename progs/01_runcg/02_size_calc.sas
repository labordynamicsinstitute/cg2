*** $Id$ ;
*** $URL$ ;
***;
*** Original author: Kevin McKinney ;
***;
*** Get information required to calculate problem size ***;
*** Calculate step1 information required to run cg ***;
*** Binary file will be output in a separate program ***;

%include './config_param.sas';

*** Variables to fill ***;

%let NN=0;
%let NCELLS=0;
%let NNPERS=0;
%let NNFIRM=0;

/*BEGINCCC
	 Get the sample size (NN) 
CCCEND*/

data strip01(keep=pik sein year);
   set MYWORK.wage_history_01;
run;

proc contents data=strip01 out=tmp;
run;

data tmp;
    set tmp; 
    if _n_ = 1 then call symput('NN',nobs);
run;

/*BEGINCCC
	Get the number of firms (NNFIRM) 
CCCEND*/

proc sort data=strip01(keep=pik sein) out=strip nodupkey;
  by sein pik;
run;

data strip firm(keep=firmseq);
  set strip;
    by sein;

  retain firmseq (0);

  if first.sein then do;
     firmseq=firmseq+1;
     output firm;
  end;
  output strip;
run;

proc contents data=firm out=tmp;
run;

data tmp;
    set tmp; 
    if _n_ = 1 then call symput('NNFIRM',nobs);
run;

*** Create the _lookup dataset         ***;
*** Create the cells for grouping prog ***;

proc sort data=strip nodupkey;      
  by pik sein; 
run; 

data MYWORK.strip02 
     MYWORK._lookup(keep=pik sein persseq firmseq) person(keep=persseq);
     merge strip01 strip;
     by pik sein;

  retain persseq (0);

  file "&cellout./cellsout";

  if first.pik then do;
     persseq=persseq+1;
     output person;
  end;

  if first.sein then do;
     /* create cellsout */
     put persseq 1-12 firmseq 14-26 ;
     output MYWORK._lookup;
  end;
  output MYWORK.strip02;
run;

proc contents data=MYWORK.strip02;
run;

/*BEGINCCC
 Get the number of persons (NNPERS) 
CCCEND*/

proc contents data=person out=tmp; 
run; 

data tmp;
    set tmp;
    if _n_ = 1 then call symput('NNPERS',nobs); 
run;

/*BEGINCCC
 Get the number of PIK SEIN cells (NCELLS) 
CCCEND*/

proc contents data=strip out=tmp;
run;

data tmp;
    set tmp; 
    if _n_ = 1 then call symput('NCELLS',nobs);
run;

/*BEGINCCC
 Pull all the information together and generate some output.
CCCEND*/

%put MYWORK.cgin:  &NN &NCELLS &NNPERS &NNFIRM;
