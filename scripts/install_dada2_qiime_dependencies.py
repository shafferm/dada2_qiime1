from rpy2 import robjects


def main():
    # load r functions
    r_install_packages = robjects.r['install.packages']
    r_source = robjects.r['source']

    # setup bioconductor
    r_source("http://bioconductor.org/biocLite.R")
    r_biocLite = robjects.r['biocLite']

    # install argparser
    r_install_packages("argparser", repos="http://cran.us.r-project.org")

    # install dada2
    r_biocLite('dada2')


if __name__ == "__main__":
    main()
