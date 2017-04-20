from rpy2 import robjects

# load r functions
r_install_packages = robjects.r['install.packages']
r_source = robjects.r['source']

# setup bioconductor
r_source("https://bioconductor.org/biocLite.R")
r_biocLite = robjects.r['biocLite']

# install argparser
install_packages("argparser", repos = "http://cran.us.r-project.org")

# install dada2
r_biocLite('dada2')
