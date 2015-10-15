# Copyright (C) 2012, Christof Buchbender
# BSD License (License.txt)
r""" Installation of Astrolyze using distutils.core.setup.
Basically it copies all scripts to
"""

import os
import sys
import site
from distutils.core import setup

r""" This function sets up the optional astrolyze databases for the maps
classes and *MAYBE LATER* also the dictionary containing the Information on
individual molecules for the ``lte`` package.
"""
from pysqlite2 import dbapi2 as sqlite
import os

def get_line_parameter(filein, database):
    r""" Reads in the line names and frequencies from ``filein`` and creates a
    table Lines in the ``database``.
    """
    const_c = 299792458.  # Speed of light [m]
    filein = open(filein).readlines()
    lines = []
    for row in filein[1:]:
        line_name, frequency = row.split()
        lines += [[line_name, float(frequency) * 1e9,
                   float(const_c / float(frequency) / 1e9)]]
    print database
    connection = sqlite.connect(database)
    cursor = connection.cursor()
    cursor.execute('CREATE TABLE Lines (id INTEGER PRIMARY KEY,'
                   'Name VARCHAR(50), '
                   'Frequency FLOAT, '
                   'Wavelenght Float)')
    for i in lines:
        cursor.execute('INSERT INTO Lines VALUES (null, ?, ?, ?)', (i[0], i[1],
                       i[2]))
    connection.commit()
    cursor.close()
    connection.close()


def get_galaxy_parameter(filein, database):
    r"""
    """
    filein = open(filein).readlines()
    galaxies = []
    for row in filein[1:]:
        (galaxy_name, morphology_type, distance, v_lsr, RA, DEC, PA,
        inclination, R25) = row.split()
        galaxies += [[galaxy_name, morphology_type, float(distance),
                     float(v_lsr), RA, DEC, float(PA), float(inclination),
                     float(R25)]]
    connection = sqlite.connect(database)
    cursor = connection.cursor()
    cursor.execute('CREATE TABLE Galaxies (id INTEGER PRIMARY KEY, Name '
                   'VARCHAR(50), MorphType VARCHAR(50), Distance DOUBLE, VLSR '
                   'DOUBLE, Central Position VARCHAR(50), PA DOUBLE, '
                   'Inclination FLOAT, R25 FLOAT)')
    for i in galaxies:
        cursor.execute('INSERT INTO Galaxies VALUES (null, ?, ?, ?, ?, ?, ?, '
                       '?, ?)',(i[0], i[1], i[2], i[3], i[4] + ' ' + i[5], i[6],
                                i[7], i[8]))
    connection.commit()
    cursor.close()
    connection.close()


def get_calibration_error(filein, database):
    r"""
    """
    filein = open(filein).readlines()
    calibration_error_list = []
    for row in filein[1:]:
        items = row.split()
        telescope = items[0]
        species = items[1]
        calibration_error = items[2]
        # The rest of the words in the row are interpreted as reference.
        # ' '.join() produces one string with one space between the items.
        reference = ' '.join(items[3:])
        calibration_error_list += [[telescope, species,
                                    float(calibration_error), reference]]
    connection = sqlite.connect(database)
    cursor = connection.cursor()
    cursor.execute('CREATE TABLE cal_error (id INTEGER PRIMARY KEY, Telescope '
                    'VARCHAR(50), Species VARCHAR(50), uncertainty DOUBLE, '
                    'Reference VARCHAR(50))')
    for i in calibration_error_list:
        cursor.execute('INSERT INTO cal_error VALUES (null, ?, ?, ?, ?)',
                       (i[0], i[1], i[2], i[3]))
    connection.commit()
    cursor.close()
    connection.close()


def create_database(database):
    r"""
    """
    if os.path.isfile(database):
        os.system('rm -rf ' + database)
    filein = os.path.expanduser('cfg/line_parameter.txt')
    get_line_parameter(filein, database)
    filein = os.path.expanduser('cfg/galaxy_parameter.txt')
    get_galaxy_parameter(filein, database)
    filein = os.path.expanduser('cfg/calibration_error.txt')
    get_calibration_error(filein, database)


create_database(os.path.expanduser('astrolyze/database/parameter.db'))


setup(
    name='astrolyze',
    version='0.1.0',
    author='Christof Buchbender',
    author_email='christof.buchbender@gmail.com',
    packages=['astrolyze',
              'astrolyze/maps',
              'astrolyze/spectra',
              'astrolyze/sed',
              'astrolyze/lte',
              'astrolyze/functions',
              'astrolyze/database'
             ],
    package_data={'astrolyze': ['database/parameter.db']},
    url='',
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
)

# change_permissions()
