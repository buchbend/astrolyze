# Copyright (C) 2012, Christof Buchbender
# BSD License (License.txt)
import math
import os
import string
import sys
import pyfits
import pywcs

import numpy as np
from pysqlite2 import dbapi2 as sqlite
from scipy.ndimage import gaussian_filter

import astrolyze_prefix as prefix
import astrolyze.functions.constants as const
from astrolyze.functions import units


class Map:
    '''
    Parent Class with functions common to all three data Formats fits, Gildas
    and Miriad.
    '''
    def __init__(self, mapName, nameConvention=True):
        '''
        Initialize a map to maps.
        '''
        self.gildas_formats = ['gdf', 'mean', 'velo', 'width', 'lmv',
                               'lmv-clean']
        self.fits_formats = ['fits']
        self.miriad_formats = ['']
        self.class_formats = ['30m', 'apex']
        self.nameConvention = nameConvention
        self.mapName = mapName
        # Test if the file exists. Directory for Miriad.
        # File for fits or GILDAS.
        if (not os.path.isdir(self.mapName)
            and not os.path.isfile(self.mapName)):
            print 'Exiting: ' + self.mapName + ' does not exist'
            sys.exit()
        # Get Informations from the Name Convention
        if nameConvention:
            self.mapNameList = mapName.split('_')
            self.comments = []
            self.source = self.mapNameList[0].split('/')[-1]
            self.prefix = self.mapNameList[0].replace(self.source, '')
            self.telescope = self.mapNameList[1]
            self.species = self.__resolveSpecies()
            self.fluxUnit = self.mapNameList[3]
            # Test for dataFormat.
            if self.mapName.endswith('.fits'):
                self.dataFormat = 'fits'
                self.mapNameList[-1] = self.mapNameList[-1].replace('.fits',
                                                                    '')
            for i in self.gildas_formats:
                if self.mapName.endswith('.' + i):
                    self.dataFormat = i
                    self.mapNameList[-1] = self.mapNameList[-1].replace('.' +
                                                                        i, '')
            for i in self.class_formats:
                if self.mapName.endswith('.' + i):
                    self.dataFormat = i
                    self.mapNameList[-1] = self.mapNameList[-1].replace('.' +
                                                                        i, '')
            if os.path.isdir(self.mapName):
                # Miriad Data Format uses directories
                self.dataFormat = ''
            self.resolution = self.__resolveResolution()
            if len(self.mapNameList) > 5:
                for i in range(len(self.mapNameList) - 6):
                    self.comments += [self.mapNameList[i + 5]]
                self.comments += [self.mapNameList[-1]]
        # Only load the File if no Name Convention is given
        if not nameConvention:
            print 'Only Files with correct naming are supported!!!'
            sys.exit()
        try:
            self.connection = sqlite.connect(str(prefix.dataBase) +
                                         'parameter.db')
            self.cursor = self.connection.cursor()
            self.cursor.execute("SELECT * FROM Galaxies WHERE Name = ?",
                            (self.source.upper(),))
            self.params = self.cursor.fetchall()[0]
            self.type = self.params[2]
            self.distance = self.params[3]
            self.vlsr = self.params[4]
            self.centralPosition = self.params[5]
            self.pa = self.params[6]
            self.inclination = self.params[7]
            self.R25 = self.params[8]
            self.cursor.close()
            self.connection.close()
        except:
            self.params = None
            self.type = None
            self.distance = None
            self.vlsr = None
            self.centralPosition = None
            self.pa = None
            self.inclination = None
            self.R25 = None
        try:
            self.connection = sqlite.connect(str(prefix.dataBase) +
                                             'parameter.db')
            self.cursor = self.connection.cursor()
            self.cursor.execute("SELECT * FROM Lines WHERE Name = ?",
                                (self.species.upper(),))
            self.params = self.cursor.fetchall()[0]
            self.frequency = self.params[2]
            self.wavelenght = self.params[3]
            self.cursor.close()
            self.connection.close()
        except:
            pass
        try:
            self.connection = sqlite.connect(str(prefix.dataBase) +
                                             'parameter.db')
            self.cursor = self.connection.cursor()
            self.cursor.execute("SELECT * FROM cal_error WHERE Telescope = "
                                " ? AND Species = ?", (self.telescope.upper(),
                                                       self.species.upper()))
            self.params = self.cursor.fetchall()[0]
            self.calibrationError = self.params[3]
            self.cursor.close()
            self.connection.close()
        except:
            self.calibrationError = np.nan
        self.get_beam_size()

    def __resolveSpecies(self):
        '''
        Gets the frequency from the map name if possible.
        '''
        species = self.mapNameList[2]
        if 'mum' in species:
            try:
                self.wavelenght = float(species.replace('mum', '')) * 1e-6
                self.frequency = 299792356 / self.wavelenght
            except:
                self.frequency = np.nan
                self.wavelenght = np.nan
        elif 'mm' in species:
            try:
                self.wavelenght = float(species.replace('mm', '')) * 1e-3
                self.frequency = 299792356 / self.wavelenght
            except:
                self.frequency = np.nan
                self.wavelenght = np.nan
        elif 'GHz' in species:
            try:
                self.frequency = float(species.replace('GHz', '')) * 1e9
                self.wavelenght = 299792356 / self.frequency
            except:
                self.frequency = np.nan
                self.wavelenght = np.nan
        else:
            self.frequency = np.nan
            self.wavelenght = np.nan
        return species

    def __resolveResolution(self):
        '''
        Read the the resolution string from the map name.
        '''
        # TODO: include handling of 'uk'
        string = self.mapNameList[4]
        # Test if there is a digit after the last point.
        # To exclude file endings like .fits.gdf.
        # In this the dataFormat would be 'gdf'.
        # but self.mapNameList[4] storing the resolution
        # Still does contain points only due to numbers.
        test = string.split('.')
        x = True
        while x:
            print string
            if 'uk' in string:
                break
            try:
                float(test[-1][:1])
                x = False
            except KeyboardInterrupt:
                sys.exit()
            except:
                string = string.replace('.' + test[-1], '')
                test = test[0:-1]
        if 'uk' in string:
            return string
        # Resolve the resolution naming scheme explained above.
        if 'x' in string and 'a' in string:
            major = float(string.split('x')[0])
            minor = float(string.split('x')[1].split('a')[0])
            pa = float(string.split('x')[1].split('a')[1])

        if 'x' in string and 'a' not in string:
            major = float(string.split('x')[0])
            minorData = string.split('x')[1]
            pa = 0.0

        if 'x' not in string and 'a' in string:
            major = float(string.split('a')[0])
            minor = float(string.split('a')[0])
            paData = string.split('a')[1]

        if 'x' not in string and 'a' not in string:
            major = float(string)
            minor = float(string)
            pa = 0
        return [major, minor, pa]

    def resolutionToString(self, resolution=None):
        if resolution is None:
            if float(self.resolution[2]) == 0.0:
                if float(self.resolution[0]) == float(self.resolution[1]):
                    string = "%1.2f" % self.resolution[0]
                if float(self.resolution[0]) != float(self.resolution[1]):
                    string = ("%1.2f" % self.resolution[0] + 'x' +
                              "%1.2f" % self.resolution[1])
            if float(self.resolution[2]) != 0.0:
                if float(self.resolution[0]) == float(self.resolution[1]):
                    string = ("%1.2f" % self.resolution[0] + 'a' +
                              "%1.1f" % self.resolution[2])
                if float(self.resolution[0]) != float(self.resolution[1]):
                    string = ("%1.2f" % self.resolution[0] + 'x' +
                              "%1.2f" % self.resolution[1] + 'a' +
                              "%1.1f" % self.resolution[2])

        if resolution is not None and type(resolution) is not str:
            if float(resolution[2]) == 0.0:
                if float(resolution[0]) == float(resolution[1]):
                    string = "%1.2f" % resolution[0]
                if float(resolution[0]) != float(resolution[1]):
                    string = ("%1.2f" % resolution[0] + 'x' +
                              "%1.2f" % resolution[1])
            if float(resolution[2]) != 0.0:
                if float(resolution[0]) == float(resolution[1]):
                    string = ("%1.2f" % resolution[0] + 'a' +
                              "%1.1f" % resolution[2])
                if float(resolution[0]) != float(resolution[1]):
                    string = ("%1.2f" % resolution[0] + 'x' +
                              "%1.2f" % resolution[1] + 'a' +
                              "%1.1f" % resolution[2])
        if type(resolution) is str:
            string = resolution
        return string

    def get_beam_size(self):
        r"""
        Calulates the Beamsize in m^2 if the distance to the source is given
        if not given the PixelSize is in sterradian.

        Notes
        -----

        The formula used is:

        .. math:

            \Omega = 1.133 * FWHM(rad)^2 \cdot (Distance(m)^2)
        """
        if self.resolution != 'uk':
            if self.distance == None:
                self.beamSize = (1.133 * const.a2r ** 2 * self.resolution[0] *
                                 self.resolution[1])
            else:
                self.beamSize = (1.133 * (self.distance * const.a2r *
                                 const.pcInM) ** 2 * self.resolution[0] *
                                 self.resolution[1])
        else:
            self.beamSize = np.nan

    def changeMapName(self, source=None, telescope=None, species=None,
                      fluxUnit=None, resolution=None, comments=None,
                      dataFormat=None, prefix=None):
        '''
        This function can be used to change the names of the maps and make a
        copy of the file to the new name and/or location.
        '''
        source = source or self.source
        telescope = telescope or self.telescope
        species = species or self.species
        fluxUnit = fluxUnit or self.fluxUnit
        resolution = resolution or self.resolution
        dataFormat = dataFormat or self.dataFormat
        prefix = prefix or self.prefix
        if comments is None:
            comments = comments or self.comments
        self.source = source or self.source
        self.telescope = telescope or self.telescope
        self.species = species or self.species
        self.fluxUnit = fluxUnit or self.fluxUnit
        self.resolution = resolution or self.resolution
        self.dataFormat = dataFormat or self.dataFormat
        self.prefix = prefix or self.prefix
        if comments is not None:
            comments = self.comments + comments
            self.comments = self.comments + comments
        if len(self.comments) == 0:
            if  str(self.mapName) != (str(prefix) + str(source) + '_' +
                                      str(telescope) + '_' + str(species) + '_'
                                      + str(fluxUnit) + '_' + str(resolution) +
                                      '.' + str(dataFormat)):
                os.system('cp ' + str(self.mapName) + ' ' +
                          str(prefix) + str(source) + '_' + str(telescope) +
                          '_' + str(species) + '_' + str(fluxUnit) + '_' +
                          self.resolutionToString(self.resolution) + '.' +
                          str(dataFormat))
                self.mapName = (str(prefix) + str(source) + '_' +
                                str(telescope) + '_' + str(species) + '_' +
                                str(fluxUnit) + '_' +
                                self.resolutionToString(self.resolution) +
                                '.' + str(dataFormat))

        if len(self.comments) != 0:
            if (str(self.mapName) != str(prefix) + str(source) + '_' +
                                     str(telescope) + '_' + str(species) +
                                     '_' + str(fluxUnit) + '_' +
                                     str(resolution) + '_' +
                                    '_'.join(self.comments) + '.' +
                                     str(dataFormat)):

                os.system('cp ' + str(self.mapName) + ' ' + str(prefix) +
                          str(source) + '_' + str(telescope) + '_' +
                          str(species) + '_' + str(fluxUnit) + '_' +
                          self.resolutionToString(self.resolution) + '_' +
                          '_'.join(self.comments) + '.' + str(dataFormat))

                self.mapName = (str(prefix) + str(source) + '_' +
                                str(telescope) + '_' + str(species) + '_' +
                                str(fluxUnit) + '_' +
                                self.resolutionToString(self.resolution) +
                                '_' + '_'.join(self.comments) +
                                '.' + str(dataFormat))

    def returnName(self, source=None, telescope=None, species=None,
                   fluxUnit=None, resolution=None, comments=None,
                   dataFormat=None, prefix=None):
        '''
        Returns the Name corresponding to the Name convention. Single keywords
        can be changed.
        '''
        source = source or self.source
        telescope = telescope or self.telescope
        species = species or self.species
        fluxUnit = fluxUnit or self.fluxUnit
        resolution = resolution or self.resolution
        dataFormat = dataFormat or self.dataFormat
        prefix = prefix or self.prefix
        if comments is None:
            comments = self.comments
        elif comments is not None:
            comments = self.comments + comments
        if len(comments) == 0:
            if dataFormat:
                return (str(prefix) + str(source) + '_' + str(telescope) +
                        '_' + str(species) + '_' + str(fluxUnit) + '_' +
                        self.resolutionToString(resolution) + '.' +
                        str(dataFormat))
            else:
                return (str(prefix) + str(source) + '_' + str(telescope) +
                        '_' + str(species) + '_' + str(fluxUnit) + '_' +
                        self.resolutionToString(resolution))
        if len(comments) != 0:
            return (str(prefix) + str(source) + '_' + str(telescope) + '_' +
                    str(species) + '_' + str(fluxUnit) + '_' +
                    self.resolutionToString(resolution) + '_' +
                    '_'.join(comments) + '.' + str(dataFormat))

    def flux_conversion(self, x=None, major=None, minor=None,
                 nu_or_lambda='nu', direction=None):
        r"""
        Calulates conversion between K.km/s and Jy/beam and vise versa.

        Parameters
        ---------
        x: float [GHz]
            Wavelenght/frequency. Defaults to the frequency of the loaded map,
            i.e. self.frequency
        major: float
            Major Axis Beam (arcsec). Default None, i.e. using self.resolution.
        minor: float
            Minor Axis Beam(arcsec). Default None, i.e. using self.resolution.
        nu_or_lambda: string
            Choose type of x: frequency = ``'nu'`` or wavelenght =
            ``'lambda'``.
        direction: string
            choose conversion direction ``'kelvin_to_jansky'``
            means Kelvin to Jansky; ``'jansky_to_kelvin'`` Jansky to Kelvin.

        See Also
        --------

        Notes
        -----

        Please note that if self.frequency and self.resolution are correctly
        set, this functions does not need any input.
        """
        if direction is not None and (direction != 'kelvin_to_jansky'
           or direction != 'jansky_to_kelvin'):
            print ('Keyword Error direction has to be kelvin_to_jansky or'
                   'jansky_to_kelvin -> Exit!')
        if self.fluxUnit in ['JyB', 'Jy/Beam'] and direction is None:
            direction = 'jansky_to_kelvin'
        if self.fluxUnit in ['Tmb', 'T', 'Kkms'] and direction is None:
            direction = 'kelvin_to_jansky'
        if (self.fluxUnit not in ['JyB', 'Jy/Beam']):
            if (self.fluxUnit not in ['Tmb', 'T', 'Kkms']):
                print ('Map is not in the right units has to be Jy/beam or '
                       'Kelvin something. -> Exit!')
                sys.exit()
            sys.exit()
        if nu_or_lambda == 'lambda':
            if direction == 'jansky_to_kelvin':
                def fcon(x, major, minor):
                    return units.jansky_to_kelvin(x, major,
                                                  minor, nu_or_lambda='lambda')
            if direction == 'kelvin_to_jansky':
                def fcon(x, major, minor):
                    return units.kelvin_to_jansky(x, major,
                                                  minor, nu_or_lambda='lambda')
        if nu_or_lambda == 'nu':
            if direction == 'jansky_to_kelvin':
                def fcon(frequency, major, minor):
                    return units.jansky_to_kelvin(x, major,
                                                  minor, nu_or_lambda='nu')
            if direction == 'kelvin_to_jansky':
                def fcon(frequency, major, minor):
                    return units.kelvin_to_jansky(x, major,
                                                  minor, nu_or_lambda='nu')
        if x is None:
            if self.frequency is not np.nan:
                x = self.frequency / 1e9
            elif  self.frequency is np.nan:
                print 'No frequency information present. Can not proceed.'
        if major is None:
            major = self.resolution[0]
        if minor is None:
            minor = self.resolution[1]
        return fcon(x, major, minor)
