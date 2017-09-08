from setuptools import setup, find_packages
from setuptools.command.install import install
from rpy2 import robjects
from glob import glob

__author__ = 'shafferm'
__version__ = '0.1.0'


class CustomInstallCommand(install):
    """Customized setuptools install command."""
    def run(self):
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
        install.run(self)

setup(
    cmdclass={'install': CustomInstallCommand},
    name="dada2_qiime1",
    version=__version__,
    install_requires=['rpy2', 'biom-format', 'numpy'],
    scripts=glob("scripts/*.py"),
    packages=find_packages(),
    description="Using DADA2 with qiime 1",
    author="Michael Shaffer",
    author_email='michael.shaffer@ucdenver.edu',
    package_data={'': ['*.r', '*.R']},
    include_package_data=True,
    url="https://github.com/shafferm/dada2_qiime1/",
    download_url="https://github.com/shafferm/dada2_qiime1/tarball/%s" % __version__
)
