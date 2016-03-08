 # Copyright (C) 2012, Christof Buchbender
# BSD License (License.txt)
import math
import os
import string
import sys
import ConfigParser

import numpy as np
from pysqlite2 import dbapi2 as sqlite
from scipy.ndimage import gaussian_filter

from generaltools import log_tools

import astrolyze.functions.constants as const
from astrolyze.functions import units


class Map(object):
    '''
    ``Map`` is the parent Class for the ``maps``-package. It contains all
    functions that are common to all supported map-formats, i.e. Fits,
    GILDAS and Miriad. This class is only supposed to be called through
    the FitsMap, GildasMap, and MiriadMap classes.

    Parameters
    ----------

    map_name : string
        The name and path of the file that is to be initialized to the maps
        package.
    '''
    def __init__(self, map_name, **kwargs):
        '''
        Initialize a map to maps.
        '''
        USER = os.getenv("USER")
        config_path = "/home/{}/.astrolyze/".format(USER)
        config_file = "astrolyze.cfg"
        self.config = ConfigParser.ConfigParser()
        self.config.read("{}{}".format(config_path, config_file))
        # Definition of the unit nomenclatures.
        self.log = log_tools.init_logger(
            directory="/home/{}/.astrolyze/".format(USER),
            name="astrolyze"
        )
        USER = os.getenv("USER")
        self.database_prefix = self.config.get("General", "database_prefix")
        self.database_prefix = self.database_prefix.format(USER)

        # Load configuration from config_file
        self.jansky_beam_names = self.config.get(
            "Units", "jansky_beam_names"
        ).split(',')
        self.jansky_pixel_names = self.config.get(
            "Units", "jansky_pixel_names"
        ).split(',')
        self.tmb_names = self.config.get(
            "Units", "tmb_names"
        ).split(',')
        self.MJy_per_sterad_names = self.config.get(
            "Units", "MJy_per_sterad_names"
        ).split(',')
        self.erg_sec_pixel_names =  self.config.get(
            "Units", "erg_sec_pixel_names"
        ).split(',')
        self.erg_sec_beam_names = self.config.get(
            "Units", "erg_sec_beam_names"
        ).split(',')
        self.erg_sec_sterad_names =  self.config.get(
            "Units", "erg_sec_sterad_names"
        ).split(',')

        self.known_units = (self.jansky_beam_names + self.jansky_pixel_names +
                            self.tmb_names + self.MJy_per_sterad_names +
                            self.erg_sec_beam_names +
                            self.erg_sec_sterad_names +
                            self.erg_sec_pixel_names)
        # Definition of possible data format endings for the different
        # programs.

        self.gildas_formats = self.config.get(
            "Formats", "gildas_formats"
        ).split(',')
        self.fits_formats = self.config.get(
            "Formats", "fits_formats"
        ).split(',')
        self.miriad_formats = ['',None]
        self.class_formats = self.config.get(
            "Formats", "class_formats"
        ).split(',')
        # name_convention is not needed anymore. Only kept for backward
        # compatibality.
        self.map_name = map_name
        # Test if the file exists. Directory for Miriad.
        # File for Fits and GILDAS.
        if ((not os.path.isdir(self.map_name)
             and not os.path.isfile(self.map_name))):
            self.log.critical(
                'Exiting: {} does not exist'.format(self.map_name)
            )
            raise SystemExit
        # Get Informations from the file name.

        self.map_nameList = map_name.split('/')[-1].split('_')
        if len(self.map_nameList) < 5:
            self.log.critical(
                "This map does not follow the Naming Convention for Astrolyze"
            )
            raise SystemExit
        self.map_nameList = map_name.split('/')[-1].split('_')
        self.prefix_list = map_name.split('/')[0:-1]
        self.comments = []
        self.source = self.map_nameList[0].split('/')[-1]
        if len(self.prefix_list) > 0:
            self.prefix = string.join(self.prefix_list, '/') + '/'
        elif len(self.prefix_list) == 0:
            self.prefix = ''
        self.telescope = self.map_nameList[1]
        self.species = self._resolveSpecies()
        self.fluxUnit = self.map_nameList[3]
        # Check dataFormat.
        if self.map_name.endswith('.fits'):
            self.dataFormat = 'fits'
            self.map_nameList[-1] = self.map_nameList[-1].replace(
                '.fits', ''
            )
        for i in self.gildas_formats:
            if self.map_name.endswith('.' + i):
                self.dataFormat = i
                self.map_nameList[-1] = self.map_nameList[-1].replace(
                    '.' + i, ''
                )
        for i in self.class_formats:
            if self.map_name.endswith('.' + i):
                self.dataFormat = i
                self.map_nameList[-1] = self.map_nameList[-1].replace(
                    '.' + i, ''
                )
        if os.path.isdir(self.map_name):
            # Miriad Data Format uses directories
            self.dataFormat = None
        self.resolution = self._resolveResolution()
        # Entries after the fifth are considered comments.
        if len(self.map_nameList) > 5:
            for i in range(len(self.map_nameList) - 6):
                self.comments += [self.map_nameList[i + 5]]
            self.comments += [self.map_nameList[-1]]
        #!!!!!  TODO: Bas implementation. Try should only contain very short
        # parts of the program otherwise errors in the program are camouflaged.
        try:
            self.connection = sqlite.connect(str(self.database_prefix) +
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
            self.connection = sqlite.connect(str(self.database_prefix) +
                                             'parameter.db')
            self.cursor = self.connection.cursor()
            self.cursor.execute("SELECT * FROM Lines WHERE Name = ?",
                                (self.species.upper(),))
            self.params = self.cursor.fetchall()[0]
            self.frequency = self.params[2]
            self.wavelength = self.params[3]
            self.cursor.close()
            self.connection.close()
        except:
            pass
        try:
            self.connection = sqlite.connect(str(self.database_prefix) +
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

    def _resolveSpecies(self):
        '''
        Gets the frequency from a database on basis of the map name if
        possible.
        '''
        species = self.map_nameList[2]
        if 'mum' in species:
            try:
                self.wavelength = float(species.replace('mum', '')) * 1e-6
                self.frequency = 299792356 / self.wavelength
            except:
                self.frequency = np.nan
                self.wavelength = np.nan
        elif 'mm' in species:
            try:
                self.wavelength = float(species.replace('mm', '')) * 1e-3
                self.frequency = 299792356 / self.wavelength
            except:
                self.frequency = np.nan
                self.wavelength = np.nan
        elif 'GHz' in species:
            try:
                self.frequency = float(species.replace('GHz', '')) * 1e9
                self.wavelength = 299792356 / self.frequency
            except:
                self.frequency = np.nan
                self.wavelength = np.nan
        else:
            self.frequency = np.nan
            self.wavelength = np.nan
        return species

    def _resolveResolution(self):
        '''
        Reads the resolution string from the map name.
        '''
        # TODO: include handling of 'uk'
        string = self.map_nameList[4]
        # Test if there is a digit after the last point.
        # To exclude file endings like .fits.gdf.
        # In this the dataFormat would be 'gdf'.
        # but self.map_nameList[4] storing the resolution
        # Still does contain points only due to numbers.
        test = string.split('.')
        x = True
        print string
        while x:
            if 'uk' in string:
                break
            try:
                float(test[-1][:1])
                x = False
            except KeyboardInterrupt:
                sys.exit()
            except:
                string = string.replace('.{}'.format(test[-1]), '')
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
        r""" Converts the resolution list to a string to be printed and
        included in the file names.
        """
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
        r""" Calculates the beam-size in steradians and in m^2. Fot the latter
        the distance to the source has to be given.

        Returns
        -------

        Initialization if the variables:
        self.beamSizeSterad and self.beamSizeM2

        Notes
        -----

        The formula used is:

        .. math:

            \Omega = 1.133 * FWHM(rad)^2 \cdot (Distance(m)^2)

        """
        if self.resolution != 'uk':
            self.beamSizeSterad = (1.133 * const.a2r ** 2 * self.resolution[0]
                                   * self.resolution[1])
            if self.distance is not None:
                self.beamSizeM2 = (1.133 * (self.distance * const.a2r *
                                            const.pcInM) ** 2 *
                                   self.resolution[0] * self.resolution[1])
        else:
            self.beamSize = np.nan

    def change_map_name(self, source=None, telescope=None, species=None,
                        fluxUnit=None, resolution=None, comments=None,
                        dataFormat=False, prefix=None):
        '''
        This function can be used to change the names of the maps and make a
        copy of the file to the new name and/or location.
        '''
        source = source or self.source
        telescope = telescope or self.telescope
        species = species or self.species
        fluxUnit = fluxUnit or self.fluxUnit
        resolution = resolution or self.resolution
        if dataFormat is False:
            dataFormat = self.dataFormat
        prefix = prefix or self.prefix
        if comments is None:
            comments = comments or self.comments

        # Now update the variables
        self.source = source
        self.telescope = telescope
        self.species = species
        self.fluxUnit = fluxUnit
        self.resolution = resolution
        if dataFormat is not False:
            self.dataFormat = dataFormat
        self.prefix = prefix 
        if comments is not None:
            comments = self.comments + comments
            self.comments = self.comments + comments
        comment_string = _comment_to_string()
        target_file_name = "{}{}_{}_{}_{}_{}{}.{}".format(
            prefix,
            source, telescope, species, fluxUnit,
            self.resolutionToString(self.resolution),
            comment_string,
            dataFormat
        )
        if len(self.comments) == 0:
            if  str(self.map_name) != (str(prefix) + str(source) + '_' +
                                       str(telescope) + '_' + str(species) +
                                       '_' + str(fluxUnit) + '_' +
                                       str(resolution) + '.' +
                                       str(dataFormat)):
                copy_command = "cp {} {}".format(self.map_name,
                                                 target_file_name)
                subprocess.call(copy_command, shell=True)
                self.map_name = target_file_name

        if len(self.comments) != 0:
            if ((str(self.map_name) != str(prefix) + str(source) + '_' +
                 str(telescope) + '_' + str(species) + '_' + str(fluxUnit) +
                 '_' + str(resolution) + '_' + '_'.join(self.comments) + '.' +
                 str(dataFormat))):
                os.system('cp ' + str(self.map_name) + ' ' + str(prefix) +
                          str(source) + '_' + str(telescope) + '_' +
                          str(species) + '_' + str(fluxUnit) + '_' +
                          self.resolutionToString(self.resolution) + '_' +
                          '_'.join(self.comments) + '.' + str(dataFormat))

                self.map_name = (str(prefix) + str(source) + '_' +
                                 str(telescope) + '_' + str(species) + '_' +
                                 str(fluxUnit) + '_' +
                                 self.resolutionToString(self.resolution) +
                                 '_' + '_'.join(self.comments) +
                                 '.' + str(dataFormat))

    def _comment_to_string(self):
        """ Converts the comment list to a string 
        """
        comment_string = "_".join(self.comments)
        if comment_string != "":
            comment_string = "_{}".format(comment_string)
        return comment_string

    def returnName(self, source=None, telescope=None, species=None,
                   fluxUnit=None, resolution=None, comments=None,
                   dataFormat=False, prefix=None):
        '''
        Returns the Name corresponding to the Name convention. Single keywords
        can be changed.

        This function is useful to generate a writeout name for a changed file
        without overwriting the current ``self.map_name``.

        Parameters
        ----------

        All possible parameters from the "Naming Convention" plus the new
        prefix.
        '''
        source = source or self.source
        telescope = telescope or self.telescope
        species = species or self.species
        fluxUnit = fluxUnit or self.fluxUnit
        resolution = resolution or self.resolution
        prefix = prefix or self.prefix
        if dataFormat is False:
            dataFormat = self.dataFormat
        if comments is None:
            comments = self.comments
        elif comments is not None:
            comments = self.comments + comments
        if len(comments) == 0:
            if dataFormat is not None:
                return (str(prefix) + str(source) + '_' + str(telescope) +
                        '_' + str(species) + '_' + str(fluxUnit) + '_' +
                        self.resolutionToString(resolution) + '.' +
                        str(dataFormat))
            if dataFormat is None:
                return (str(prefix) + str(source) + '_' + str(telescope) +
                        '_' + str(species) + '_' + str(fluxUnit) + '_' +
                        self.resolutionToString(resolution))
        if len(comments) != 0:
            if dataFormat is not None:
                return (str(prefix) + str(source) + '_' + str(telescope) +
                        '_' + str(species) + '_' + str(fluxUnit) + '_' +
                        self.resolutionToString(resolution) + '_' +
                        '_'.join(comments) + '.' + str(dataFormat))
            if dataFormat is None:
                return (str(prefix) + str(source) + '_' + str(telescope) +
                        '_' + str(species) + '_' + str(fluxUnit) + '_' +
                        self.resolutionToString(resolution) + '_' +
                        '_'.join(comments))

    def flux_conversion(self, x=None, major=None, minor=None,
                        nu_or_lambda='nu', direction=None):
        r"""
        Calulates conversion between K.km/s and Jy/beam and vise versa.

        Parameters
        ----------

        x : float [GHz]
            Wavelength/frequency. Defaults to the frequency of the loaded map,
            i.e. self.frequency
        major : float
            Major Axis Beam (arcsec). Default None, i.e. using self.resolution.
        minor : float
            Minor Axis Beam(arcsec). Default None, i.e. using self.resolution.
        nu_or_lambda : string
            Choose type of x: frequency = ``'nu'`` or wavelength =
            ``'lambda'``.
        direction : string
            choose conversion direction ``'kelvin_to_jansky'``
            means Kelvin to Jansky; ``'jansky_to_kelvin'`` Jansky to Kelvin.

        Notes
        -----

        If self.frequency and self.resolution are correctly set, this functions
        does not need any input. Otherwise this has to be given explicitly.
        """
        if direction is not None and (direction != 'kelvin_to_jansky'
                                      or direction != 'jansky_to_kelvin'):
            print ('Keyword Error direction has to be kelvin_to_jansky or'
                   'jansky_to_kelvin -> Exit!')
        if self.fluxUnit in ['JyB', 'Jy/Beam'] and direction is None:
            direction = 'jansky_to_kelvin'
        if self.fluxUnit in ['Tmb', 'T', 'Kkms'] and direction is None:
            direction = 'kelvin_to_jansky'
        print self.fluxUnit
        if ((self.fluxUnit not in ['JyB', 'Jy/Beam'] and self.fluxUnit
             not in ['Tmb', 'T', 'Kkms'])):
                print ('Map is not in the right units has to be Jy/beam or '
                       'Kelvin something. -> Exit!')
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
