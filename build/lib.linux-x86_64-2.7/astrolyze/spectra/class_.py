# Copyright (C) 2012, Christof Buchbender
# BSD Licencse
import math
import os
import gc
import string
import sys
import pyfits
import pywcs
import subprocess

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
    def __init__(self, map_name, nameConvention=True):
        r"""Initializes a Class Spectral File."""
        astrolyze.maps.main.Map.__init__(self, map_name, nameConvention)
        self._init_file_to_class()
        if self.dataFormat not in self.class_formats:
            print ('Exiting: Not a CLASS format (AFAIK). Supported'
                   'Formats Can be extended.')
            sys.exit()

    def _init_file_to_class(self):
        r"""
        Initial setup, loading all the variables correponding to
        the cube.
        """
        self.set_defaults()
        pyclass.comm('file in ' + self.map_name)
        pyclass.comm('find')
        pyclass.comm('get first')
        self._load_class_variables()
        self.dec_coordinate = self.vars_.r.head.pos.bet.__sicdata__
        self.ra_coordinate = self.vars_.r.head.pos.lam.__sicdata__
        # conversion to degrees
        self.ra_coordinate = self.ra_coordinate * const.r2d
        self.dec_coordinate = self.dec_coordinate * const.r2d
        self.central_coordinate_degrees = [self.ra_coordinate,
                                           self.dec_coordinate]
        self.central_coordinate_equatorial = astFunc.degrees_to_equatorial(
                                             self.central_coordinate_degrees)
        self.found = self.vars_.found.__sicdata__
        pyclass.comm('find')

    def _load_class_variables(self):
        r"""
        This function actually imports the Class variables.
        """
        self.vars_ = pyclass.gdict

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
        print "Central Coordinate", self.central_coordinate_equatorial
        print "calculated offset", offset
        pyclass.comm('file in ' + self.map_name)
        pyclass.comm('set offset ' + str(offset[0]) + ' ' + str(offset[1]))
        while True:
            pyclass.comm('set match ' + str(accuracy))
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

    def get_circular_region_from_cube(self, coordinates, aperture,
                                      angle=0, prefix=None):
        r""" Read a all spectra falling into a circle with diameter aperture
        and center coordinate.

        Parameters
        ----------

        coordinate: list
            Equatorial coordinates, e.g. [['1:34:7.00', '+30:47:52.00']]

        aperture: float
            The diameter of the circle in arcsec.

        angle: float
            If the cube was rotated before the angle has to be specified
            to calculate the correct offset.

        prefix: string
            The new path where the averaged spectrum will be stored.

        Returns
        -------
        A new 30m file with only the spectra inside the circle.
        """

        if prefix is None:
            _prefix = self.prefix
        if prefix is not None:
            _prefix = prefix
        offsets = []
        pyclass.comm('file in ' + self.map_name)
        self._init_file_to_class()
        for coord in coordinates:
            print coord
            coords = coord[1]
            source = coord[0]
            print source, coords
            offsets += [[source, [astFunc.calc_offset(self.central_coordinate_equatorial,
                                         coords, angle, output_unit='arcsec')]]]
        return_name = self.returnName(prefix = _prefix,
                                      comments=['all' ,'circle'])
        pyclass.comm('las\\file out "' + return_name + '" single /overwrite' )
        print '#######', offsets
        pyclass.comm('find')
        self._load_class_variables()
        found = self.vars_.found.__sicdata__
        for i in range(found):
            pyclass.comm("get n")
            self._load_class_variables()
            bet_off_0 = float(str(self.vars_.off_beta.__sicdata__))
            lambd_off_0 = float(str(self.vars_.off_lambda.__sicdata__))
            for offset in offsets:
                bet_off_new = bet_off_0 - offset[1][0][1]
                lambd_off_new = lambd_off_0 - offset[1][0][0]
                if ((math.floor(math.sqrt((bet_off_new) ** 2 +
                                      (lambd_off_new) ** 2)) <
                     (aperture/2))):
                    pyclass.comm("mod source " + str(offset[0]))
                    pyclass.comm("write")
                gc.collect()
            gc.collect()
        pyclass.comm("file in " + return_name)
        old_name = return_name
        for offset in offsets:
            return_name = self.returnName(prefix = _prefix,
                                      comments=[str(offset[0]) ,'circle'])
            pyclass.comm('las\\file out "' + return_name + '" single /overwrite' )
            pyclass.comm('find /source ' + str(offset[0]))
            pyclass.comm('copy')
        subprocess.call('rm ' + old_name, shell=True)

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

    def go_where(self, name=False):
        r"""
        Helper function that saves the current plot.
        """
        name = name or self.returnName(dataFormat='eps')
        fileout = open("dummy.class","w")
        fileout.write('file in ' + self.map_name + '\n'
                      'find\n'
                      'go where\n'
                      r'GTVL\ha ' + name + ' /dev eps color /over'+'\n'
                      'exit\n')
        fileout.close()
        os.system('class -nw @dummy.class')
        subprocess.call('rm dummy.class',shell=True)


    def save_figure(self, name=None):
        r"""
        Helper function that saves the current plot.
        """
        name = name or self.returnName(dataFormat='eps')
        print r'GTVL\ha ' + name + ' /dev eps color /over'
        pyclass.comm(r'GTVL\ha ' + name + ' /dev eps color /over')

    def quick_view(self, number=1, window=True):
        r"""
        Helper Functions that displays the first spectrum of the loaded
        file.
        """
        # TODO: make more usefull.
        pyclass.comm('file in ' + self.map_name)
        if window:
            pyclass.comm('dev im w')
        pyclass.comm('find')
        pyclass.comm('get {}'.format(number))
        pyclass.comm('pl')

    def get_freq_at_max(self, range_=None):
        r''' Get the frequency at the maximum of the spectrum
        '''
        if range_:
            self.max_rx = list(self.ry).index(np.max(self.ry[range_[0]:range_[1]]))
        if not range_:
            self.max_rx = list(self.ry).index(np.max(self.ry))
        self.max_freq = self.freqs[self.max_rx]

    def get_freq_array(self, dic):
        r''' Calculates an array of frequencies corresponding to the individual
        channels.

        !!!! --- NOT FINISHED --- !!!!
        '''
        self.freqs = []
        for i in range(dic.channels):
            self.freqs += [self.central_frequency - ((self.refchan-(i+1)) *
                                                     self.freq_step)]
        return self.freqs

    def get(self, number='next'):
        '''Wrapper to the get function with the possibility to get additional
        information on the spectrum.
        '''
        if ((type(number) == str) and 
            (number in ['next', 'first', 'last', 'n' ,'f' ,'l'])):
            pyclass.comm('get {}'.format(number))
        elif type(number) == int:
            pyclass.comm('get {}'.format(number))
        else:
            raise TypeError
        self.rx = pyclass.gdict.rx.__sicdata__
        self.ry = pyclass.gdict.ry.__sicdata__
        self.refchan = pyclass.gdict.reference.__sicdata__
        self.central_frequency = pyclass.gdict.frequency
        self.freq_step = pyclass.gdict.freq_step
        self.freq_array = self.get_freq_array(pyclass.gdict)
        self.get_freq_at_max()
    def make_noise_table(self, sigma_tolerance=5,
                         lower_window='-5', upper_window='5'):
        pyclass.comm('file in ' + self.map_name)
        pyclass.comm('find')
        pyclass.comm('set u v')
        data_vars = pyclass.gdict
        found = data_vars.found.__sicdata__
        pyclass.comm('set window ' +  lower_window + ' ' + upper_window)
        model_template = []
        for i in range(found):
            i = i + 1
            pyclass.comm('get n')
            pyclass.comm('bas 0')
            data_vars = pyclass.gdict
            sigma = data_vars.sigma.__sicdata__
            off_beta = float(str(data_vars.off_beta.__sicdata__))
            off_lambda = float(str(data_vars.off_lambda.__sicdata__))
            x_pc = self.distance * const.a2r * off_beta
            y_pc = self.distance * const.a2r * off_lambda
            channels = data_vars.channels.__sicdata__
            print sigma
            good_channels = []
            for i in range(channels):
                i  = i + 1
                if i > 1:
                    channel_before = data_vars.ry[i-1].__sicdata__
                else:
                    channel_before = None
                if i < channels:
                    channel_after = data_vars.ry[i+1].__sicdata__
                else:
                    channel_after = None
                ry_i = data_vars.ry[i].__sicdata__
                if ry_i >= 5*sigma:
                    if (channel_before >= sigma_tolerance*sigma or
                        channel_after >= sigma_tolerance*sigma):
                        channel = data_vars.channel.__sicdata__
                        velocity = data_vars.velocity.__sicdata__
                        velo_step = data_vars.velo_step.__sicdata__
                        chan_vel = velocity+((channel-i)*velo_step)
                        good_channels += [chan_vel]
            if good_channels != []:
                for chan in good_channels:
                    model_template += [(x_pc, y_pc, chan)]
            good_channels = []
            print '\n\n'
        return model_template

    def fit_hfs_hcn(self, line_window=[None,None], number=None, line=None,
                    figure_path='tmp/'):
        r''' Work in Progress.....
        '''
        # Create the HCN input parameters
        fileout = open('hfs-hcn.dat', 'w')
        fileout.write('3\n')
        fileout.write('-7.064 0.2\n')
        fileout.write('0  1\n')
        fileout.write('4.842  0.6\n')
        fileout.close()
        # Fit has to be done for the velocity axis
        pyclass.comm('set u v')
        # Window for base-lining
        pyclass.comm('set window {} {}'.format(line_window[0],line_window[1]))
        pyclass.comm('find')
        data_vars = pyclass.gdict
        if not number:
            for i in range(1, self.found,1):
                pyclass.comm('get n')
                pyclass.comm('bas 0')
                pyclass.comm('method hfs hfs-hcn.dat')
                try:
                    pyclass.comm('min')
                except:
                    continue
                pyclass.comm('pl')
                pyclass.comm('vis')
                line_name = str(data_vars.line.__sicdata__)
                pyclass.comm('ha ' + figure_path + '_' + line_name +
                             str(i)  +  '.eps /dev eps color /over')
               # intensities = data_vars.r.head.hfs.result[0] * relative intensities
                #positions = data_vars.r.head.hfs.result[1] +- relative positions
                fwhm = data_vars.r.head.hfs.result[2]
                print fwhm
        if number:
            pyclass.comm('get {}'.format(number))
            pyclass.comm('bas 0')
            pyclass.comm('method hfs hfs-hcn.dat')
            pyclass.comm('min')
            pyclass.comm('pl')
            pyclass.comm('vis')
        subprocess.call('rm hfs-hcn.dat',shell=True)


    def fit_gauss(self, line_window=[None,None], number=None, lines=None,
                  figure_path='tmp/'):
        r''' Work in Progress.....
        '''
        if not os.path.isdir(figure_path):
            subprocess.call('mkdir ' + figure_path, shell=True)
        pyclass.comm('set u v')
        pyclass.comm('find')
        if not number:
            results = []
            data_vars = pyclass.gdict
            for i in range(1, self.found,1):
                line_string = ('lines  ' + str(len(lines)) +  ' "')
                for j in lines:
                    line_string += '" "'.join(j) + '" '
                line_string += '/nocursor'
                print line_string
                pyclass.comm('get n')
                pyclass.comm('set window {} {}'.format(line_window[0],line_window[1]))
                pyclass.comm('bas 0')
                pyclass.comm('method gauss')
                pyclass.comm(line_string)
                try:
                    pyclass.comm('min')
                except:
                    continue
                line_name = str(data_vars.line.__sicdata__)
                pyclass.comm('pl')
                pyclass.comm('vis')
                pyclass.comm('ha ' + figure_path + line_name + '_' + str(i) +
                             '.eps /dev eps color /over')
                intensities = [float(data_vars.r.head.gau.nfit[0].__sicdata__),
                               float(data_vars.r.head.gau.nerr[0].__sicdata__)]
                positions = [float(data_vars.r.head.gau.nfit[1].__sicdata__),
                             float(data_vars.r.head.gau.nerr[1].__sicdata__)]
                fwhm = [float(data_vars.r.head.gau.nfit[2].__sicdata__),
                        float(data_vars.r.head.gau.nerr[2].__sicdata__)]
                results += [[i, intensities, positions, fwhm]]
        if number:
            pyclass.comm('get {}'.format(number))
            pyclass.comm('bas 0')
            pyclass.comm('method gauss')
            pyclass.comm('min')
            pyclass.comm('pl')
            pyclass.comm('vis')
            intensities = [data_vars.r.head.gau.nfit[0],
                           data_vars.r.head.gau.nerr[0]]
            positions = [data_vars.r.head.gau.nfit[1],
                         data_vars.r.head.gau.nerr[1]]
            fwhm = [data_vars.r.head.gau.nfit[2], data_vars.r.head.gau.nerr[2]]
            results = [number, intensities, positions, fwhm]
        return results





    # def toFits(self, freq, mol_name, vel_range=295, overwrite=True,
    #               no_extract=False):
    #     pyclass.comm("file in " + self.map_name))
    #     if os.path.isfile(fileout):
    #         pyclass.comm("file out {}".format(self.returnName().replace('30m')))
    #     if not os.path.isfile(fileout) or overwrite:
    #         pyclass.comm("file out {} single /over".format(fileout))
    #     pyclass.comm("find")
    #     vars = pyclass.gdict
    #     found = vars.found.__sicdata__
    #     for i in range(found):
    #         pyclass.comm("get n")
    #         pyclass.comm("mod freq {}".format(freq))
    #         source = str(vars.source.__sicdata__).strip()
    #         if not no_extract:
    #             pyclass.comm("extract 'velocity-{0}' "
    #                          "'velocity+{0}'".format(vel_range/2))
    #         subprocess.call('rm ../fits_spectra/M33{}-{}.'
    #                         'fits'.format(source, mol_name),
    #                         shell=True)
    #         print ("fits write ../fits_spectra/M33{}-{}"
    #                " /mode spec".format(source, mol_name))
    #         pyclass.comm("fits write ../fits_spectra/M33{}-{}"
    #                " /mode spec".format(source, mol_name))
    #         pyclass.comm("write")
