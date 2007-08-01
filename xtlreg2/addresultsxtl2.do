cap mata: mata drop addresultscgxtl2()
mata

function addresultscgxtl2(string scalar dependent, string rowvector covariates,
				string scalar individualid, string scalar unitid, string scalar pastunitid,
				real scalar nind, real scalar nunit, 
				string scalar filerawoutput,
				string scalar prefix) {


	real scalar ncov;
	real scalar fp;
	real vector betas, indeffect, uniteffect;
	real vector params;
	real matrix data;
	real scalar n,i;
	string scalar indeffect_name, uniteffect_name, pastuniteffect_name;
	string scalar cmd_listcovnames;
	string scalar cmd_listcovvalues;

	ncov = length(covariates);

	printf("addresultscgxtl : adding results to dataset\n");
	printf("dependent: %s\n", dependent);
	printf("individualid: %s\n", individualid);
	printf("unitid: %s\n",unitid);
	printf("pastunitid: %s\n", pastunitid);
	printf("nind: %g\n", nind);
	printf("nunit: %g\n", nunit);
	printf("filerawoutput: %s\n",filerawoutput);
	printf("prefix: %s\n", prefix);
	printf("ncov: %g\n", ncov);
	
	fp = fopen(filerawoutput,"r");
	params = fgetmatrix(fp);
	fclose(fp);
	
	betas 		= params[1 		      	.. ncov				];
	indeffect 		= params[ncov + 1 		.. ncov + nind			];
	uniteffect 		= params[ncov + nind + 1 	.. ncov + nind + nunit -1	];

	indeffect_name 		= "indeffect";
	uniteffect_name 	= "uniteffect";
	pastuniteffect_name 	= "pastuniteffect";

	st_addvar("double", indeffect_name);
	st_addvar("double", uniteffect_name);
	st_addvar("double", pastuniteffect_name);	

	st_view(data, . ,(individualid, unitid, pastunitid, indeffect_name, uniteffect_name, pastuniteffect_name));

	n = rows(data);

	for (i = 1; i<= n ; i ++) {
		data[i,4] = indeffect[data[i,1]];
		if (data[i,2] != nunit) {
			data[i,5] = uniteffect[data[i,2]];
		} else {
			data[i,5] = 0;
		}
		if (data[i,3] != nunit && data[i,3] != .) {
			data[i,6] = uniteffect[data[i,3]];
		} else if (data[i,3] == nunit) {
			data[i,6] = 0;
		} else {
			data[i,6] = .;
		}
	}

	stata("gen xb = 0");

	cmd_listcovnames  =	"";
	cmd_listcovvalues = 	"";

	for (i = 1; i<= ncov; i++) {
		cmd_listcovnames 	= sprintf("%s %s",cmd_listcovnames, covariates[i]);
		cmd_listcovvalues	= sprintf("%s %g,",cmd_listcovvalues, betas[i]);
		stata(sprintf("replace xb = xb + %g * %s",betas[i], covariates[i]));
	}

	stata(sprintf("matrix input betas = (%s)",cmd_listcovvalues));
	stata(sprintf("matrix colnames betas = %s",cmd_listcovnames));
	
}

end
