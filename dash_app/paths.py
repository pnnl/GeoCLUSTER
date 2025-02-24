import os
# from dotenv import load_dotenv
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

localpath = ""
internal_pnnl_path = "/var/www/html/dash_app/"
aws_path = "/www/GeoCLUSTER/dash_app/"

# load_dotenv()
print(os.environ)
deployment_type = os.getenv("deployment_type", "aws")
absolute_path = localpath if deployment_type == "local" else aws_path

inpath_dict = {

		"absolute_path": absolute_path,
	
	
	# A local URL prefix for file requests. 
	# Defaults to url_base_pathname, and must end with routes_pathname_prefix. 
	# env: DASH_REQUESTS_PATHNAME_PREFIX.

	# Author's Recommendation: If get a screen that just says "Loading..." then need to 
	# 						 manipulate the following.
	

		"requests_pathname_prefix": "/", # "/GeoCLUSTER/"


	
	# A local URL prefix to use app-wide. Default '/'. 
	# Both requests_pathname_prefix and routes_pathname_prefix default to url_base_pathname. 
	# env: DASH_URL_BASE_PATHNAME

	# Author's Recommendation: Not needed and will often throw an error but consider editing 
	# 						 if requests_pathname_prefix alone is not working.
	

		"url_base_pathname": "/", # "/GeoCLUSTER/"


	
	# Data inpaths and/or outpaths.
	

	"h5_filepath": absolute_path + "data/clgs_results_final_float32.h5", # cuts memory of the data in half

	"geoCLUSTER_results_pathname": absolute_path + 'tmp/geoCLUSTER_results.xlsx',

	"properties_H2O_pathname": absolute_path + "data/properties_H2O.mat",
	"properties_CO2v2_pathname": absolute_path + "data/properties_CO2v2.mat",
	"additional_properties_CO2v2_pathname": absolute_path + "data/additional_properties_CO2v2.mat",

	"tmatrix_pathname": absolute_path + 'data/tmatrix.csv',

	"summary_tbl_in": absolute_path + "data/summary_tbl_in.csv"

}





