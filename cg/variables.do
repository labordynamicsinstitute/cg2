program defvar

	mata:	veccov = J(1,1,("keystage1","keystage2"))
	mata: veclabel = J(1,1,("Key Stage 1 Dummy","Key Stage 2 Dummy"))	

	global DEPENDENT = "y"
	global COVARIATES = "keystage1 keystage2"
	global NCOV = 1

end
