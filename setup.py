#!/usr/bin/env python

from setuptools import setup

import os


def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'refinem', 'VERSION'))
    return versionFile.readline().strip()

if __name__ == '__main__':

    dirName = os.path.dirname(__file__)
    if dirName and os.getcwd() != dirName:
        os.chdir(dirName)

    #f = open('README.rst', 'rt')
    #long_description = f.read()
    #f.close()

    setup(
        name='refinem',
        version=version(),
        author='Donovan Parks',
        author_email='donovan.parks@gmail.com',
        packages=['refinem', 'refinem.plots'],
        scripts=['bin/refinem'],
        package_data={'refinem' : ['VERSION', './distributions/*.txt','./data_files/hmms/*.hmm']},
        url='http://pypi.python.org/pypi/refinem/',
        license='GPL3',
        description='A toolbox for improving population genomes.',
        #long_description=long_description,
        install_requires=[
            "numpy>=1.9.0",
            "matplotlib>=1.4.0",
            "biolib>=0.0.45",
            "jinja2>=2.7.3",
            "mpld3>=0.2",
	    "pysam",
            "weightedstats"],
    )
