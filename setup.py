# Copyright (C) 2012, Christof Buchbender
# BSD License (License.txt)
r""" Installation of Astrolyze using distutils.core.setup.
Basically it copies all scripts to
"""

import os
import sys
import site
from distutils.core import setup
setup(
    name='astrolyze',
    version='0.1.2',
    author='Christof Buchbender',
    author_email='buchbend@ph1.uni-koeln.de',
    url='https://github.com/buchbend/astrolyze.git',
    packages=['astrolyze',
              'astrolyze/maps',
              'astrolyze/spectra',
              'astrolyze/sed',
              'astrolyze/lte',
              'astrolyze/functions',
              'astrolyze/database'
             ],
    data_files = [("cfg/",
                   ["cfg/calibration_error.txt",
                    "cfg/galaxy_parameter.txt",
                    "cfg/line_parameter.txt"])],
    license='LICENSE.txt',
    description=('Reduction and analysing tools for (mainly)'
                 'Radioastronomical Data.'),
    long_description=open('README.txt').read(),
    classifiers=['Topic :: Scientific/Engineering :: Astronomy'],
    requires=[
        "numpy",
        "pyfits",
        "matplotlib",
        "scipy",
        "pywcs",
        "pysqlite2",
        "docutils"
    ],
    scripts=['scripts/astrolyze_opt_db_setup.py']
)

# change_permissions()
