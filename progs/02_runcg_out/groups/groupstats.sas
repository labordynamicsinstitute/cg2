*** get some simple stats on the groups ***;


libname dot '.';

data groups;
   infile 'groups';

   length persseq firmseq group 5;
   input persseq firmseq group;
run;

proc freq data=groups order=freq;
   tables group /out=dot.groupfreq;
run;

proc freq data=dot.groupfreq;
   table count;
run;

proc print data=dot.groupfreq(obs=100);
run;
