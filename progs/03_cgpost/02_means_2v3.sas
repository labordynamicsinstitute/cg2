libname temp '.';
 
data temp.means02;
infile '02_means_2v3.dat' lrecl=28;
input
learn
exper
;
 
proc contents;
run;
 
proc print;
run;
