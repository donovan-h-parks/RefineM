from distutils.core import setup

import os


def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'refinem', 'VERSION'))
    return versionFile.read().strip()

setup(
    name='refinem',
    version=version(),
    author='Donovan Parks',
    author_email='donovan.parks@gmail.com',
    packages=['refinem', 'refinem.plots'],
    scripts=['bin/refinem'],
    package_data={'refinem' : ['VERSION', './distributions/*.txt']},
    url='http://pypi.python.org/pypi/refinem/',
    license='GPL3',
    description='A toolbox for improving population genomes.',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy >= 1.9.0",
        "matplotlib >= 1.4.0",
        "biolib >= 0.0.19",
        "jinja2 >= 2.7.3"
        "mpld3 >= 0.2"],
)
