#!/usr/bin/env python

# Copyright (C) 2009-2012, Christof Buchbender
# BSD Licencse
r""" This function sets up the optional astrolyze databases for the maps
classes and *MAYBE* also the dictionary containing the Information on
individual molecules for the ``lte`` package.
"""
from pysqlite2 import dbapi2 as sqlite
import os
import astrolyze.functions.constants as const


def get_line_parameter(filein, database):
    r""" Reads in the line names and frequencies from ``filein`` and creates a
    table Lines in the ``database``.
    """
    filein = open(filein).readlines()
    lines = []
    for row in filein[1:]:
        line_name, frequency = row.split()
        lines += [[line_name, float(frequency) * 1e9,
                   float(const.c / float(frequency) / 1e9)]]
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
    filein = os.path.expanduser('~/.astrolyze/setup/line_parameter.txt')
    get_line_parameter(filein, database)
    filein = os.path.expanduser('~/.astrolyze/setup/galaxy_parameter.txt')
    get_galaxy_parameter(filein, database)
    filein = os.path.expanduser('~/.astrolyze/setup/calibration_error.txt')
    get_calibration_error(filein, database)


if __name__ == "__main__":
    print os.path.expanduser('~/.astrolyze/database/parameter.db')
    os.system('ls ' + os.path.expanduser('~/.astrolyze/database/parameter.db'))
    create_database(os.path.expanduser('~/.astrolyze/database/parameter.db'))
