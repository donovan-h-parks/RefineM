from distutils.core import setup

import os

def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'cleanm', 'VERSION'))
    return versionFile.read().strip()

setup(
    name='cleanm',
    version=version(),
    author='Donovan Parks',
    author_email='donovan.parks@gmail.com',
    packages=['cleanm'],
    scripts=['bin/cleanm'],
    package_data={'cleanm' : ['VERSION']},
    url='http://pypi.python.org/pypi/comparem/',
    license='GPL3',
    description='A toolbox for improving population genomes.',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy >= 1.8.0",
        "biolib >= 0.0.1"],
)
