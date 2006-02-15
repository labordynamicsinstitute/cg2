*** Create some fake PIK SEIN records to use as test data for CG2 ***; 

options obs=max ls=180 ps=10000 mprint mlogic symbolgen;

%include "./config_param.sas";

data cg.wage_history_01;

  length pik $2 sein $3 year 4;

  input pik $ sein $ year  learn exper;


  datalines;
12 555 1991 10.045 22
12 555 1992 10.095 23
12 555 1993 10.135 24
12 555 1994 10.139 24.5
12 625 1995 10.215 25
12 625 1996 10.235 26
12 625 1997 10.225 27
12 625 1998 10.335 28
12 625 1999 10.325 29
16 455 1996 10.805 41
16 455 1997 10.845 42
16 455 1998 10.965 42.75
16 455 1999 11.045 43.5
16 725 1991 11.145 44
16 725 1992 11.235 45
16 725 1993 11.245 46
16 725 1994 11.245 47
16 725 1995 11.255 48
20 555 1992 11.345 45
20 555 1993 11.465 46
20 555 1994 11.545 47
20 555 1995 11.675 48
20 825 1996 10.945 49
20 825 1997 11.145 50
20 825 1998 11.225 51
20 825 1999 11.235 52
24 555 1994 10.805 34
24 555 1995 10.915 34.5
24 555 1996 10.925 35.2
24 555 1997 10.945 36
28 625 1996 10.545 29
28 625 1997 10.575 30
28 625 1998 10.625 31
28 625 1999 10.665 32
32 725 1992 11.445 52
32 725 1993 11.465 53
32 725 1994 11.475 54
32 725 1995 11.465 55
36 455 1996 10.020 24
36 455 1997 10.035 25
36 455 1998 10.040 26
36 455 1999 10.045 27
40 825 1991 10.805 42
40 825 1992 10.935 43
40 825 1993 10.175 44
40 825 1994 11.435 45
40 825 1995 10.955 46
40 825 1996 11.165 47
40 825 1997 11.295 48
;
run;

proc contents;
proc means;
proc print;
run;
