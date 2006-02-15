*** Create a file similar to cellsout ***;
*** This file contains the cells in cellsout sorted by firm ***;

libname dot '.';

libname pcells '../';

options obs=max;

proc sort data=pcells._lookup(keep=firmseq persseq) out=temp0;
    by firmseq persseq;
run;

data _NULL_;
   set temp0;

   file 'firmcells';

   put persseq 1-12 firmseq 14-26 ;
run;
