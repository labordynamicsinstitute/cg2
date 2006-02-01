*************************************************************************************
*														*
*    This file outputs a nice LaTeX table with the results provided in CGOUT		*
*														*
************************************************************************************* 


cd "W:\projets\schooleffects\src\cg"
do params.do

** First step : merge effects with the original cgin data file

mata

//*** Creating dataset with pupil effects

fp = fopen("cgout_$STANDARDIZE_$TRANSFORMATION_$DEPENDENT","r")
params = fgetmatrix(fp);
pupileffects = params[$NCOV+1 .. $NCOV + $NPUPILS];
st_addvar("double","pupileffect")
st_addvar("float","pupilid")
stata("set obs $NPUPILS")
stata("replace pupilid = _n")
st_store((1,$NPUPILS),"pupileffect",pupileffects)
stata("sort pupilid")
stata("save pupileffects.dta, replace")
fclose(fp)
stata("clear")

//*** Creating dataset with school effects

fp=fopen("cgout_$STANDARDIZE_$TRANSFORMATION_$DEPENDENT","r")
params = fgetmatrix(fp);
schooleffects = params[$NCOV+$NPUPILS+1 .. $NCOV + $NPUPILS + $NSCHOOLS -1];
st_addvar("double","schooleffect")
st_addvar("float","schoolid")
neffects = $NSCHOOLS-1
stata("set obs "+strofreal(neffects))
stata("replace schoolid = _n")
st_store((1,$NSCHOOLS-1),"schooleffect",schooleffects)
stata("sort schoolid")
stata("save schooleffects.dta, replace")
fclose(fp)
stata("clear")

end

use cgin_englishmaths, clear
sort pupilid
merge pupilid using pupileffects
tab _merge

drop _merge
sort schoolid
merge schoolid using schooleffects
tab _merge
replace schooleffect = 0 if  _merge == 1 

//*** Adding covariates : label & value

defvar

mata
fp = fopen("cgout_$STANDARDIZE_$TRANSFORMATION_$DEPENDENT","r")
params = fgetmatrix(fp);
fclose(fp)

betas = params[1 .. $NCOV];

for (i = 1; i<= $NCOV; i++) {
	stata(sprintf("gen beta_%s = %f", veccov[i], params[i]));
	stata(sprintf("label variable %s %s %s %s ", veccov[i], char(34),veclabel[i],char(34)));
}
end

gen xb = 0

foreach covvar in $COVARIATES {
	replace xb = xb + beta_`covvar' * `covvar'
}

save cgout_final, replace

// *** Creating LaTeX file with table (means, variance and correlations)

defvar

mata

st_view(y=.,.,"$DEPENDENT")
st_view(cov=., 	., veccov)
st_view(pupileffect=., 	., "pupileffect")
st_view(schooleffect=., ., "schooleffect")
st_view(variableseffect=., ., "xb")

pupileff_mean 	= mean(pupileffect,1)
pupileff_stddev 	= variance(pupileffect,1)

schooleff_mean 	= mean(schooleffect,1)
schooleff_stddev 	= variance(schooleffect,1)

y_mean 		= mean(y,1)
y_stddev 		= variance(y,1)

variableseff_mean = mean(variableseffect,1)
variableseff_stddev = variance(variableseffect,1)

corr_y_pupileff   = ((1/length(y))*y'pupileffect - y_mean*pupileff_mean)/sqrt(pupileff_stddev*y_stddev)
corr_y_schooleff   = ((1/length(y))*y'schooleffect - y_mean*schooleff_mean)/sqrt(schooleff_stddev*y_stddev)
corr_pupileff_schooleff = ((1/length(y))*schooleffect'pupileffect - schooleff_mean*pupileff_mean)/sqrt(pupileff_stddev*schooleff_stddev)
corr_vars_pupileff = ((1/length(y))*variableseffect'pupileffect - variableseff_mean*pupileff_mean)/sqrt(pupileff_stddev*variableseff_stddev)
corr_vars_schooleff = ((1/length(y))*variableseffect'schooleffect - variableseff_mean*schooleff_mean)/sqrt(schooleff_stddev*variableseff_stddev)
corr_vars_y = ((1/length(y))*variableseffect'y - variableseff_mean*y_mean)/sqrt(y_stddev*variableseff_stddev)

fp=fopen("W:\projets\schooleffects\tables\decvar_$STANDARDIZE_$TRANSFORMATION_$DEPENDENT.tex","w")

fwrite(fp,sprintf("\\begin{table}\r\n"));
fwrite(fp,sprintf("\\caption{Summary Statistics for the Decomposition of Variance for Individual Data, Maths,  Key Stage~1 1998-1999 and Key Stage~2 2002-2003}\r\n"));
fwrite(fp,sprintf("\\begin{center}\r\n \\begin{tabular}{lcccccc}\r\n \hline \hline\r\n"));
fwrite(fp,sprintf("& & & \\multicolumn{4}{c}{Simple Correlation with:} \\\ \r\n"));
fwrite(fp,sprintf("\\multicolumn{1}{c}{Variable Description} & Mean & Std.~Dev. & \$y\$ & \$x \\beta\$ & \$\\theta\$ & \$\\psi\$ \\\ \r\n"));
fwrite(fp,sprintf("\hline \\\ \r\n"));
fwrite(fp,sprintf("\$y\$, Standardized of Maths Grade & %9.4f & %9.4f & 1.0000 & %9.4f & %9.4f & %9.4f \\\ \r\n",y_mean,y_stddev,corr_vars_y,corr_y_pupileff,corr_y_schooleff));
fwrite(fp,sprintf("\$x\\beta\$, Predicted Effect of \$x\$ variables & %9.4f & %9.4f & %9.4f & 1.0000 & %9.4f & %9.4f \\\ \r\n", variableseff_mean, variableseff_stddev, corr_vars_y,corr_vars_pupileff,corr_vars_schooleff));
fwrite(fp,sprintf("\$\\theta\$, Pupil Effect & %9.4f & %9.4f & %9.4f & %9.4f & 1.0000 & %9.4f  \\\ \r\n",pupileff_mean,pupileff_stddev,corr_y_pupileff,corr_vars_pupileff,corr_pupileff_schooleff));
fwrite(fp,sprintf("\$\\psi\$, School Effect & %9.4f & %9.4f &  %9.4f & %9.4f & %9.4f & 1.0000 \\\ \r\n",schooleff_mean,pupileff_stddev,corr_y_schooleff,corr_vars_schooleff,corr_pupileff_schooleff));
fwrite(fp,sprintf("\\\ \r\n"));
fwrite(fp,sprintf("\\hline \r\n"));
fwrite(fp,sprintf("\\end{tabular} \r\n"));
fwrite(fp,sprintf("\\end{center} \r\n"));
fwrite(fp,sprintf("\\sourcenpd \r\n"));
fwrite(fp,sprintf("\\end{table} \r\n"));

fclose(fp)

end
