library(reticulate)
library(rgee)

# Get the username
HOME <- Sys.getenv("HOME")

# Install miniconda
reticulate::install_miniconda()

# Set global parameters
Sys.setenv("RETICULATE_PYTHON" = sprintf("%s/Library/r-miniconda/bin/python3", HOME))
reticulate::py_config()
reticulate::py_discover_config()
py_install("numpy")

# install rgee Python dependencies
ee_install(py_env = "rgee")
# retart R session 

library(rgee)
ee_Initialize() 
# this will open up a webpage to generate a token for notebook access. 
# use wenjing.xuuu@gmail.com 

# update EE crediential 
ee_Authenticate(auth_mode = "notebook")

ee_check() # checked 
ee_Initialize(). #maybe I just need to log in every time. 
