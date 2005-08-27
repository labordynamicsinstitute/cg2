libname temp '.';
 
data temp.rhs;
infile '01_rhs.dat' lrecl=14;
input
exper
;
 
proc contents;
run;
 
proc print;
run;
