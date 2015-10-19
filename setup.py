# Copyright (C) 2012, Christof Buchbender
# BSD License (License.txt)
r""" Installation of Astrolyze using distutils.core.setup.
Basically it copies all scripts to
"""

import os
import sys
import site
from distutils.core import setup


SUDO_USER = os.getenv("SUDO_USER")
print SUDO_USER
setup(
    name='astrolyze',
    version='0.1.3',
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
    data_files = [("/home/{}/.astrolyze/".format(SUDO_USER),
                   ["cfg/calibration_error.txt",
                    "cfg/galaxy_parameter.txt",
                    "cfg/line_parameter.txt"])],
    license='LICENSE.txt',
    description=('Reduction and analysing tools for (mainly)'
                 'Radioastronomical Data.'),
    long_description=open('README.txt').read(),
    requires=[
        "numpy",
        "pyfits",
        "matplotlib",
        "scipy",
        "pywcs",
        "pysqlite2",
        "docutils",
        "generaltools"
    ],
    classifiers=[
          'Intended Audience :: Science/Research',
          'License :: BSD License',
          'Operating System :: OS Independent',
          'Programming Language :: C',
          'Programming Language :: Python :: 2.7',
          'Topic :: Scientific/Engineering :: Astronomy',
          'Topic :: Scientific/Engineering :: Physics'
      ],
    scripts=['scripts/astrolyze_opt_db_setup.py']
)

# change_permissions()
