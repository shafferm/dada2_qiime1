from setuptools import setup, find_packages
from setuptools.command.install import install
from rpy2 import robjects

__author__ = 'shafferm'


class CustomInstallCommand(install):
    """Customized setuptools install command - prints a friendly greeting."""
    def run(self):
        # load r functions
        r_install_packages = robjects.r['install.packages']
        r_source = robjects.r['source']

        # setup bioconductor
        r_source("https://bioconductor.org/biocLite.R")
        r_biocLite = robjects.r['biocLite']

        # install argparser
        r_install_packages("argparser", repos="http://cran.us.r-project.org")

        # install dada2
        r_biocLite('dada2')
        install.run(self)

setup(
    cmdclass={'install': CustomInstallCommand},
    name="dada2_qiime1",
    version="0.0.0",
    scripts=["scripts/qiime_dada2.py"],
    packages=find_packages(),
    description="Using DADA2 with qiime 1",
    author="Michael Shaffer",
    author_email='michael.shaffer@ucdenver.edu',
    package_data={'': ['*.r', '*.R']},
    include_package_data=True,
    install_requires=['rpy2', 'biom-format', 'numpy']
)
