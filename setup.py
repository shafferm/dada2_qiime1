from setuptools import setup, find_packages
from glob import glob

__author__ = 'shafferm'
__version__ = '0.1.2'


setup(
    name="dada2_qiime1",
    version=__version__,
    install_requires=['rpy2 ==2.8.5', 'biom-format', 'numpy', 'qiime'],
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
