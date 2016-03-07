# Copyright (C) 2012, Christof Buchbender
# BSD Licencse
import math
import os
import string
import sys

from pysqlite2 import dbapi2 as sqlite
from scipy.ndimage import gaussian_filter
import pyclass

from astrolyze.maps.main import *
import astrolyze.maps.gildas

import astrolyze.functions.constants as const
from astrolyze.functions import astro_functions as astFunc
import astrolyze.functions.units

class ClassSpectra(Map):
    r"""
    Provides some usefull automated functions to work on Class
    Spectra in an convenient way.

    Examples
    --------

    Extracting a spectra at a given position from a spectral cube can be done
    as follows

    >>> from astrolyze.spectra import *
    >>>
    >>> cube = ClassSpectra(filename)
    >>> coordinate = ['1:34:7.00', '+30:47:52.00']
    >>> cube.get_spectra_from_cube(coordinate)
    Generates a 30m file with comment extract in the actual cube.prefix path.
    """
    def __init__(self, map_name):
        r"""Initializes a Class Spectral File."""
        astrolyze.maps.main.Map.__init__(self, map_name)
        print self.map_name
        self._init_file_to_greg()
        if self.dataFormat not in self.class_formats:
            print ('Exiting: Not a CLASS format (AFAIK). Supported'
                   'Formats Can be extended.')
            sys.exit()

    def _init_file_to_greg(self):
        r"""
        Initial setup, loading all the variables correponding to
        the cube.
        """
        self.set_defaults()
        pyclass.comm('file in ' + self.map_name)
        pyclass.comm('find')
        pyclass.comm('get first')
        self._load_class_variables()
        self.dec_coordinate = self.vars.r.head.pos.bet.__sicdata__
        self.ra_coordinate = self.vars.r.head.pos.lam.__sicdata__
        # conversion to degrees
        self.ra_coordinate = self.ra_coordinate * const.r2d
        self.dec_coordinate = self.dec_coordinate * const.r2d
        self.central_coordinate_degrees = [self.ra_coordinate,
                                           self.dec_coordinate]
        self.central_coordinate_equatorial = astFunc.degrees_to_equatorial(
                                             self.central_coordinate_degrees)

    def _load_class_variables(self):
        r"""
        This function actually imports the Class variables.
        """
        self.vars = pyclass.gdict

    def set_defaults(self):
        r"""
        Reset all selection criteria.
        """
        pyclass.comm('set def')
        pyclass.comm('clear')
        pyclass.comm('set plot hist')
        pyclass.comm('set form long')
        pyclass.comm('set angle sec')

    def get_spectra_from_cube(self, coordinate, angle=0, prefix=None,
                              accuracy=2, region=False):
        r"""
        Extracts one spectra at the position of coordinates from a
        spectral cube.

        Parameters
        ----------

        coordinate: list
            Equatorial coordinates, e.g. ['1:34:7.00', '+30:47:52.00']

        angle: float
            If the cube was rotated before the angle has to be specified
            to calculate the correct offset.

        prefix: string
            The new path where the averaged spectrum will be stored.

        accuracy: float
            The tolerance in arcsec to find a spectra corresponding to the
            given coordinate.

        region: True or False
            Returns either all spectra found ``True`` or only the first
            ``False``.

        Returns
        -------

        30m file
            With the first spectrum in the list of spectra within the accuracy
            range with the given coordinate.
        """
        if prefix is None:
            _prefix = self.prefix
        if prefix is not None:
            _prefix = prefix
        offset = astFunc.calc_offset(self.central_coordinate_equatorial,
                                     coordinate, angle, output_unit='arcsec')
        self.set_defaults()
        print self.central_coordinate_equatorial
        print offset
        pyclass.comm('file in ' + self.map_name)
        self.log('file in ' + self.map_name)
        pyclass.comm('set offset ' + str(offset[0]) + ' ' + str(offset[1]))
        self.log('set offset ' + str(offset[0]) + ' ' + str(offset[1]))
        while True:
            pyclass.comm('set match ' + str(accuracy))
            self.log('set match ' + str(accuracy))
            pyclass.comm('find')
            try:
                pyclass.comm('get f')
            except KeyboardInterrupt:
                raise KeyboardInterrupt
            except:
                print  '###\nNothing found, raising tolerance by 1 arsec.\n###'
                accuracy = accuracy + 1
            else:
                print ('###\nFound a spectra in a ' + str(accuracy) + ' arcsec '
                       'radius.\n###')
                break
        if not region:
            print 'asdsa',_prefix
            return_name = self.returnName(prefix = _prefix,
                                          comments=['extract'])
            pyclass.comm('file out ' + return_name + ' single /overwrite')
            pyclass.comm('write')
        if region:
            return_name = self.returnName(prefix = _prefix, comments=['region'])
            pyclass.comm('file out ' + return_name + ' single /overwrite')
            pyclass.comm('find')
            pyclass.comm('copy')
        return ClassSpectra(return_name)

#    def set_selection(self, telescope=None, line=None, source=None):
#        r"""
#        Select subsets of the spectra in the input file.
#        """
#        telescope = telescop

    def get_region_from_cube(self, coordinate, angle=0, prefix=None, 
                             accuracy=10):
        r"""
        The same as :py:func:``get_spectra_from_cube`` but returns all spectra
        found inside a circular region arounf coordinate and in a radius of 
        accuracy arcsec. ("set match "'accuracy')
        """
        return self.get_spectra_from_cube(coordinate, angle=angle,
                                          prefix=prefix, accuracy=accuracy,
                                          region=True)

    def get_intensity_in_velocity_range(self, vel_range):
        pyclass.comm("get f")
        central_velocity = pyclass.gdict.reference.__sicdata__
        reference_channel = pyclass.gdict.reference.__sicdata__
        velocity_channel_increment = pyclass.gdict.velo_step.__sicdata__
        velocity_channel_increment = pyclass.gdict.velo_step.__sicdata__
        pyclass.comm("")
        return intensity

    def get_average_spectrum(self, prefix=None):
        r"""
        Averages all spectra in a cube.

        Parameters
        ----------

        prefix: string
            The new path where the averaged spectrum will be stored.

        Notes
        -----

        So far no selection is made so the files of the input file have to be
        consistent.
        """
        _prefix =  prefix or self.prefix
        self.set_defaults()
        pyclass.comm('file in ' + self.map_name)
        pyclass.comm('find')
        pyclass.comm('get f')
        pyclass.comm('set no match')
        pyclass.comm('aver')
        return_name = self.returnName(prefix = _prefix, comments=['average'])
        pyclass.comm('file out ' + return_name + ' single /overwrite')
        pyclass.comm('write')
        return ClassSpectra(return_name)

    def save_figure(self, name=None):
        r"""
        Helper function that saves the current plot.
        """
        name = name or self.returnName(dataFormat='eps')
        pyclass.comm('ha ' + name + '/dev eps color')

    def quick_view(self, number=1):
        r"""
        Helper Functions that displays the first spectrum of the loaded
        file. 
        """
        # TODO: make more usefull.
        pyclass.comm('file in ' + self.map_name)
        pyclass.comm('dev im w')
        pyclass.comm('find')
        pyclass.comm('get {}'.format(number))
        pyclass.comm('pl')
