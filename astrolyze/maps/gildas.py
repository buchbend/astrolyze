# Copyright (C) 2012, Christof Buchbender
# BSD Licencse
import math
import os
import string
import sys
import pyfits
import pywcs

from pysqlite2 import dbapi2 as sqlite
from scipy.ndimage import gaussian_filter
import pygreg

import main
import fits
import miriad

from astrolyze.spectra import *
from astrolyze.setup.paths import prefix
import astrolyze.functions.constants as const


class GildasMap(main.Map):
    r"""
    Wrapping GILDAS functions to use them inline with Python.
    """
    def __init__(self, mapName, nameConvention=True):
        r"""Initializes a Gildas map"""
        main.Map.__init__(self, mapName, nameConvention)
        self.__init_map_to_greg()
        if self.dataFormat not in self.gildas_formats:
            print ('Exiting: Not a Gildas format (AFAIK). Supported'
                   'Formats Can be extended.')
            sys.exit()

    def __init_map_to_greg(self):
        pygreg.comm('image ' + self.mapName)
        self._load_gildas_variables()

    def _load_gildas_variables(self):
        r"""
        This function load the Gildas variables and
        extracts from some of the SicVar instances variables consiting of
        sicname, sicdata and siclevel, to just the data part to be used by
        the python parts of the programs. The extracted variables are to be
        extended while adding functionality to the GildasMap mehtods.
        """
        self.vars = pygreg.gdict
        # load the map dimensions stored in g_dim 
        self.dimensions = self.vars.g_dim.__sicdata__
        self.naxis_1 = self.dimensions[0]
        self.naxis_2 = self.dimensions[1]
        # Load the array with, central Pixel, value, and increment
        self.convert = self.vars.g_convert.__sicdata__
        self.crpix_1 = self.convert[0][0]
        self.crval_1 = self.convert[0][1]
        self.cdelt_1 = self.convert[0][2]
        self.crpix_2 = self.convert[1][0]
        self.crval_2 = self.convert[1][1]
        self.cdelt_2 = self.convert[1][2]

    def regrid_to_arcsec(self, value):
        r"""
        Regrids the pixel size of the map to a multiple of arcseconds.

        Parameters
        ----------
        value: float
            The new pixel size in arcsecs.

        Notes
        -----
        .. warning::
            Old function no guarantee of functionality. Test or remove!
        """
#        fitsFile = self.toFits()
#        self.naxis1 = float(fitsFile.header['naxis1'])
#        self.naxis2 = float(fitsFile.header['naxis1'])
#        self.cdelt1 = float(fitsFile.header['CDELT1'])
#        self.cdelt2 = float(fitsFile.header['CDELT2'])
#        self.crval1 = float(fitsFile.header['CRVAL1'])
#        self.crval2 = float(fitsFile.header['CRVAL2'])
#        self.crpix1 = float(fitsFile.header['CRPIX1'])
#        self.crpix2 = float(fitsFile.header['CRPIX2'])
        self.cdelt1arcs = float(self.cdelt_11) / (value * (1. / 60 / 60))
        self.naxis1New = self.naxis_1 * math.sqrt(self.cdelt1arcs ** 2)
        self.crpix1New = self.naxis1New / (self.naxis_1 / self.crpix_1)
        self.comments += ['IP1']
        __interpolate_init = open('interpolate.init', 'w')
        __init_string = ('TASK\FILE \"Input file\" Y_NAME$ ' +
                     str(self.mapName) + '\n'
                     'TASK\FILE \"Output file\" X_NAME$ ' +
                     str(self.returnName()) + '\n'
                     'TASK\INTEGER \"Number of pixels\" NX$ ' +
                     str(int(self.naxis1New)) + '\n'
                     'TASK\REAL \"New reference pixel\" REFERENCE$ ' +
                     str(int(self.crpix1New)) + '\n'
                     'TASK\REAL \"New value at reference pixel\" VALUE$ ' +
                     str(self.crval_1) + '\n'
                     'TASK\REAL \"New increment\" INCREMENT$ ' +
                     str(self.cdelt1arcs) + ' \n'
                     'TASK\GO\n')
        __interpolate_init.write(__init_string)
        __fileout.close()
        __fileout = open('interpolate.greg', 'w')
        __fileout.write('run interpolate interpolate.init /nowindow\nexit\n')
        __fileout.close()
        os.system('greg -nw @interpolate.greg')

    def spectrum(self, coordinates, fileout=None, prefix=None,):
        r"""
        Wrapper to the spectrum command from greg that extracts
        a spectrum from a cube at a given positions.

        Parameters
        ----------

        coordinates: list
            A list with the coordinates in floats in units of Degrees,
            or in string for equatorial coordinates.

        fileout: string
            The name of the table where the spectrum will be stored.
            Default is the same name as the map with ``".tab"`` as ending.

        prefix:
            The path to the folder where the newly created file will be
            stored.  Defaults to the prefix currently stored in self.prefix.

        Notes
        -----

        Tested and working.

        Examples
        --------

        >>> from astrolyze.maps import *
        >>> map =  GildasMap('M33_PdBI_12co10_Tmb_22.0_2kms.gdf')
        >>> coordinate = ['1:34:7.00', '+30:47:52.00']
        >>> map.spectrum(coordinate)
        Creates M33_PdBI_12co10_Tmb_22.0_2kms.tab in the present folder.
        """
        os.system('rm spectrum.init')
        if prefix is None:
            prefix = self.prefix
        if fileout is None:
            fileout = self.returnName(prefix=prefix, dataFormat='tab')
        else:
            fileout = prefix + fileout
        __maskInit = open('spectrum.init', 'w')
        __init_string = ('! Spectrum.init\n'
                      'TASK\FILE "Input file name" IN$ "' + self.mapName +
                      '"\nTASK\FILE "Output table name" OUT$ "'
                      + fileout
                      + '"\n'
                      'TASK\CHARACTER "Absolute coordinates" COORD$ "'
                      + coordinates[0] + ' ' + coordinates[1] + '"\n'
                      'TASK\GO')
        __maskInit.write(__init_string)
        __maskInit.close()
        __convFile = open('temp.greg', 'w')
        __convFile.write('run spectrum spectrum.init /nowindow\nexit\n')
        __convFile.close()
        os.system('greg -nw @temp.greg')
        os.system('rm temp.greg')
        os.system('rm spectrum.init')

    def lmv(self, fileout=None, prefix=None):
        r"""
        Wrapper to the lmv command of Class to extract spectra from a
        spectral cube.

        Parameters
        ----------

        fileout: string
            The name of the class file to write the spectra to. Defaults to
            the mapname with .30m ending.

        prefix: string
            The path were the class file will be stores. Defaults to
            the current path.

        Notes
        -----

        Tested and working.
        """
        if prefix is None:
            prefix = self.prefix
        if fileout is None:
            fileout = self.returnName(prefix=prefix, dataFormat='30m')
        else:
            fileout = prefix + fileout
        string = ('file out ' + str(fileout) + ' single /overwrite\n'
                  'lmv ' + str(self.mapName) + '\n'
                  'exit')
        classScript = open('temp.class', 'w')
        classScript.write(string)
        classScript.close()
        os.system('class -nw @temp.class')
        os.system('rm temp.class')
        return ClassSpectra(fileout)

    def mask(self, polygon, prefix=None):
        r"""
        Wrapper to the GREG task mask:

        Parameters
        ----------

        polygon: string
            path to a GILDAS polygon file with ending ``".pol"``

        prefix: string
            The path where the output is to be stored if different
            from the current prefix of the map.

        Returns
        -------

        mapObject: The masked map object.

        Examples
        --------

        >>> map.mask('poly/sourceA.pol')
        """
        if prefix is None:
            prefix = self.prefix
        __maskInit = open('mask.init', 'w')
        __comment = [polygon.split('/')[-1].replace('.pol', '')]
        __init_string = ('TASK\FILE "Polygon data file" POLYGON$ "' +
                     str(polygon) + '" \n'
                     'TASK\FILE "Input image" Y_NAME$ "' +
                     str(self.mapName) + '"\n'
                     'TASK\FILE "Output image" X_NAME$ "'
                     + str(self.returnName(prefix=prefix,
                     comments=__comment)) + '"\n'
                     'TASK\LOGICAL "Mask inside (.true.) or outside '
                     '(.false.)" MASK_IN$ NO\n'
                     'TASK\LOGICAL "Modify blanking value" MODIFY$ NO\n'
                     'SAY "If YES, Fill information Below"\n'
                     'TASK\REAL "New blanking value" BLANKING$ 0 \n'
                     'TASK\GO')
        __maskInit.write(__init_string)
        __maskInit.close()
        os.system('more mask.init')
        __convFile = open('temp.greg', 'w')
        __convFile.write('run mask mask.init /nowindow\nexit\n')
        __convFile.close()
        os.system('greg -nw @temp.greg')
        os.system('rm temp.greg')
        os.system('rm mask.init')
        return GildasMap(self.returnName(prefix=prefix,  comments=__comment))

    def reproject(self, template=None, coord=None, prefix=None,
                  keep_pixsize=False):
        r"""
        Wraps the GREG task reproject. Either use template *or* coord.

        Parameters
        ----------

        template: string
            Full path to a map in GDF Format whose central
            coordinate and pixel size will serve as a template.

        coord: list
            List of coordinate strings in RA DEC (J2000) that
            will become the new centre of the map.

        prefix: string
            The path where the output is to be stored if different
            from the current prefix of the map. If None the current
            self.prefix of the GildasMap instance is used.

        keep_pixsize: bool
            If False reproject guesses the new pixel sizes after reprojection
            these are normally smaller than the original ones.
            If True the old pixel sizes are enforced.

        Returns
        -------

        GildasMap Object: Instance for the reprojected map.

        Raises
        ------

        SystemExit
            If both **template** and **coord** are not ``None``.

        ValueError
            If keep_pixsize is not a boolean.

        Examples
        --------

        >>> map.reproject(coord = ['1:34:32.8', '30:47:00.6'])
        >>> map.reproject(template = 'M33_SPIRE_250_JyB_18.1.gdf')

        References
        ----------
        For more information on the Gildas task see:
        .. [1] www.iram.fr/GILDAS/

        """
        
        if template is not None and coord is not None:
            print "Please use either template OR coord."
            raise SystemExit
        if type(keep_pixsize) != bool:
            raise ValueError('keep_pixsize has to be True or False.')
        __reproInit = open('reproject.init', 'w')
        if prefix is None:
            prefix = self.prefix
        if template:
            comment = ['repro']
            __init_string = ('TASK\FILE "Input file" Z_NAME$ "' +
                          str(self.mapName) + '"\n'
                         'TASK\FILE "Output file" X_NAME$ "' +
                          str(self.returnName(prefix=prefix, comments=comment))
                         + '"\n'
                         'TASK\LOGICAL "Reproject on a third file\'s '
                         'projection [YES|NO]" TEMPLATE$ YES\n'
                         'TASK\FILE "This Template File Name" Y_NAME$ "' +
                         str(template) + '"\n'
                         'TASK\CHARACTER "Projection type" PROJECTION$ '
                         '"RADIO" /CHOICE     NONE GNOMONIC ORTHOGRAPHIC  '
                         'AZIMUTHAL STEREOGRAPHIC LAMBERT AITOFF RADIO\n'
                         'TASK\CHARACTER "Coord.System [EQUATORIAL [epoch]|'
                         'GALACTIC]" SYSTEM$ "EQUATORIAL 2000"\n'
                         'TASK\CHARACTER "Coord. 1 of projection center. '
                         'Could be UNCHANGED"     CENTER_1$ "01:34:11.8" '
                         '/CHOICE UNCHANGED *\n'
                         'TASK\CHARACTER "Coord. 2 of projection '
                         'center" CENTER_2$ "+30:50:23.4" /CHOICE UNCHANGED '
                         '*\n'
                         'TASK\REAL "Position angle of projection" ANGLE$ '
                         '0.000000000000000 /RANGE 0 360\n'
                         'TASK\INTEGER "Dimensions of output image '
                         '[0 0 mean guess]" DIMENSIONS$[2]  0 0\n'
                         'TASK\REAL "First axis conversion formula [0 0 0 '
                         'mean guess]" AXIS_1$[3]  0 0 0\nTASK\REAL "Second '
                         'axis conversion formula [0 0 0 mean guess]" '
                         'AXIS_2$[3] 0 0 0\n'
                         'TASK\LOGICAL "Change blanking value" CHANGE$ NO\n'
                         'TASK\REAL "New blanking value and tolerance" '
                         'BLANKING$[2]  0 0\nTASK\GO\n')
            __reproInit.write(__init_string)
        if coord:
            comment = ['repro']
            __init_string = ('TASK\FILE "Input file" Z_NAME$ "' +
                         str(self.mapName) + '"\n'
                         'TASK\FILE "Output file" X_NAME$ "' +
                         str(self.returnName(prefix=prefix, comments=comment))
                         + '"\n'
                         'TASK\LOGICAL "Reproject on a third file\'s '
                         'projection [YES|NO]" TEMPLATE$ NO\n'
                         'TASK\FILE "This Template File Name" Y_NAME$ ""\n'
                         'TASK\CHARACTER "Projection type" PROJECTION$ '
                         '"RADIO" /CHOICE     NONE GNOMONIC ORTHOGRAPHIC  '
                         'AZIMUTHAL STEREOGRAPHIC LAMBERT AITOFF RADIO\n'
                         'TASK\CHARACTER "Coord.System [EQUATORIAL [epoch]|'
                         'GALACTIC]" SYSTEM$ "EQUATORIAL 2000"\n'
                         'TASK\CHARACTER "Coord. 1 of projection center. '
                         'Could be UNCHANGED"     CENTER_1$ "' +
                          str(coord[0]) + '" /CHOICE UNCHANGED *\n'
                          'TASK\CHARACTER "Coord. 2 of projection center" '
                          'CENTER_2$ "' +
                          str(coord[1]) + '" /CHOICE UNCHANGED *\n'
                          'TASK\REAL "Position angle of projection" ANGLE$ '
                          '0.000000000000000 /RANGE 0 360\n')
            if keep_pixsize:
                __init_string += ('TASK\INTEGER "Dimensions of output image [0'
                                  '0 mean guess]" DIMENSIONS$[2] ' +
                                  str(self.naxis_1) + ' ' + str(self.naxis_1) +
                                  '\n' 'TASK\REAL "First axis conversion'
                                  'formula [0 0 0 ' 'mean guess]" AXIS_1$[3] ' 
                                  + str(self.crpix_1) + ' ' 
                                  + str(self.crval_1) +  ' ' 
                                  + str(self.cdelt_1) + '\n' 'TASK\REAL'
                                  '"Second axis conversion formula [0 0 0 '
                                  'mean guess]" AXIS_2$[3] ' +
                                  str(self.crpix_2) + ' ' +
                                  str(self.crval_2) + ' ' + 
                                  str(self.cdelt_2) + '\n')
            if not keep_pixsize:
                __init_string += ('TASK\INTEGER "Dimensions of output image [0'
                                  '0 mean guess]" DIMENSIONS$[2] 0 0 \n'
                                  'TASK\REAL "First axis conversion formula '
                                  '[0 0 0 mean guess]" AXIS_1$[3] 0 0 0 \n '
                                  'TASK\REAL "Second axis conversion formula '
                                  '[0 0 0 mean guess]" AXIS_2$[3] 0 0 0 \n')
            __init_string += ('TASK\LOGICAL "Change blanking value" CHANGE$ '
                              'YES\n TASK\REAL "New blanking value and '
                              'tolerance" BLANKING$[2]  0 0\n'
                              'TASK\GO\n')
            __reproInit.write(__init_string)
        __reproInit.close()
        os.system('more temp.greg')
        # write the Greg Script
        __convFile = open('temp.greg', 'w')
        __convFile.write('run reproject reproject.init /nowindow\nexit\n')
        __convFile.close()
        # execute the Greg Script and then delete it
        os.system('greg -nw @temp.greg')
        #os.system('rm temp.greg')
        os.system('rm reproject.init')
        return GildasMap(self.returnName(prefix=prefix, comments=comment))

    def moments(self, velo_range=[0, 0], threshold=0,
               smooth='YES', prefix=None, comment=None):
        r"""
        Wraps the GREG task moments creating the first three moments
        of the map.

        Parameters
        ----------
        velo_range: list
            Velocity range for the integration.

        threshold: float
            Value under which pixels are blanked.

        smooth: string
            One of Either ``"NO"`` or ``"YES"``. Controls
            if the map is smoothed before applying the cut threshold. Getting
            rid of noise peaks over the threshold.
            Defaults to ``'YES'``

        prefix: string
            The path where the output is to be stored if different
            from the current prefix of the map.

        comment: string
            Optional comments to be added to the new map name.

        Returns
        -------
        mean : MapObject
            The zeroth moment, i.e. the integrated intensity, is returned as a
            GildasMap object.
        """
        print self.mapName
        if comment is None:
            comment = []
        save_comments = self.comments
        self.comments = []
        if prefix is None:
            prefix = self.prefix
        __momentInit = open('moments.init', 'w')
        __init_string = ('TASK\FILE "Input file name" IN$ "' +
                       str(self.mapName) + '"\n'
                       'TASK\FILE "Output files name (no extension)" OUT$ "' +
                       self.returnName(resolution='dummy', prefix=prefix,
                           comments=[]).replace('.' + self.dataFormat,
                                                     '') + '"\n'
                       'TASK\REAL "Velocity range" VELOCITY$[2]  ' +
                        str(velo_range[0]) + ' ' + str(velo_range[1]) + '\n'
                       'TASK\REAL "Detection threshold" THRESHOLD$ ' +
                        str(threshold) + ' \n'
                       'TASK\LOGICAL "Smooth before detetction ?" SMOOTH$ ' +
                        str(smooth) + ' \n'
                       'TASK\GO\n')
        __momentInit.write(__init_string)
        __momentInit.close()
        __convFile = open('temp.greg', 'w')
        __convFile.write('run moments moments.init /nowindow\nexit\n')
        __convFile.close()
        os.system('greg -nw @temp.greg')
        os.system('mv ' + self.returnName(resolution='dummy',
                  prefix=prefix, comments=[], dataFormat='mean') + ' ' +
                  self.returnName(prefix=prefix, comments=save_comments,
                  dataFormat='mean'))
        os.system('mv ' + self.returnName(resolution='dummy',
                  prefix=prefix, comments=[], dataFormat='velo') + ' ' +
                  self.returnName(prefix=prefix, comments=save_comments,
                  dataFormat='velo'))
        os.system('mv ' + self.returnName(resolution='dummy',
                  prefix=prefix, comments=[], dataFormat='width') + ' ' +
                  self.returnName(prefix=prefix, comments=save_comments,
                  dataFormat='width'))
        #os.system('rm moments.init')
        print self.returnName(prefix=prefix,
                                         comments=save_comments,
                                         dataFormat='mean')
        return GildasMap(self.returnName(prefix=prefix,
                                         comments=save_comments,
                                         dataFormat='mean'))

    def goRot(self, angle, prefix=None):
        r"""
        Wrapper to the GREG go rot command, which rotates maps around their
        central coordinate stored in the header.

        Parameters
        ----------
        angle: float [deg]
            Rotation angle.

        prefix: string
            The path where the output is to be stored if different
            from the current prefix of the map.

        Returns
        -------
        GildasMap Object: Instance for the reprojected map.

        Examples
        --------

        >>> map.goRot(45)

        To change the central coordinate first use
        :py:func:`maps.gildas.GildasMap.reproject` e.g.:

        >>> map = map.reproject(coord=['new_RA_string','new_DEC_string'])
        >>> map.goRot(45)
        """
        __rotate = open('rotateTemp.greg', 'w')
        __rotate.write('let name ' + addOn.myStrip(str(self.mapName),
                          len(self.dataFormat) + 1) + '\n'
                          'let type ' + str(self.dataFormat) + '\n'
                          'let angle ' + str(angle) + '\n'
                          'go rot\n'
                          'exit\n')
        __rotate.close()
        os.system('greg -nw @rotateTemp.greg')
        if len(str(angle).split('.')) > 1:
            comments += [('rot' + str(str(angle).split('.')[0]) +
                               'p' + str(str(angle).split('.')[1]) + 'deg')]
        else:
            comments += ['rot' + str(angle) + 'deg']

        os.system('mv ' + addOn.myStrip(str(self.mapName),
                  len(self.dataFormat) + 1) + '-rot' + str(angle) + 'deg.'
                  + str(self.dataFormat) + ' '
                  + str(self.returnName(prefix=prefix, comments=comment)))
        os.system('rm -rf rotateTemp.greg')
        return GildasMap(self.returnName(prefix=prefix, comments=comment))

    def smooth(self, new_resolution, old_resolution=None, prefix=None):
        r"""
        Wrapper to the GREG task gauss_smooth.

        Parameters
        ----------
        new_resolution: float or list
            The resulting resolution after the smoothing.
            It can be:
                1. a float: i.e. the final major and minor beamsize.
                   The position angle will default to 0.
                2. a list with two floats: [major_axis, minor_axis]. The
                   position angle defaults to 0.
                3. a list with three floats: [major_axis, minor_axis,
                   position_angle].

        old_resolution: float or list
            Same format as new_resolution. Defaults to self.resolution of the
            map instance.

        prefix: string
            The path where the output is to be stored if different
            from the current prefix of the map.

        Notes
        -----

        .. warning::

            The gauss_smooth Task from GILDAS only gives correct output
            units when the map is on a temperature or \"per pixel\" scale.
            **Maps in Jy/Beam won't be in Jy/Beam after smoothing.**

        """
        if prefix is None:
            prefix = self.prefix
        if old_resolution == None:
            old_major = self.resolution[0]
            old_minor = self.resolution[1]
            pa = self.resolution[2]
        if old_resolution != None:
            if type(old_resolution) is list:
                if len(old_resolution) == 2:
                    old_major = old_resolution[0]
                    old_minor = old_resolution[1]
                    pa = 0
                if len(old_resolution) == 3:
                    old_major = old_resolution[0]
                    old_minor = old_resolution[1]
                    pa = old_resolution[2]
            if type(old_resolution) is not list:
                old_major = old_resolution
                old_minor = old_resolution
                pa = 0
        if type(new_resolution) is list:
            if len(new_resolution) == 2:
                new_major = new_resolution[0]
                new_minor = new_resolution[1]
                pa = 0
            if len(new_resolution) == 3:
                new_major = new_resolution[0]
                new_minor = new_resolution[1]
                pa = new_resolution[2]
        if type(new_resolution) is not list:
            new_major = new_resolution
            new_minor = new_resolution
            pa = 0
        new_resolution = [new_major, new_minor, pa]

        if (float(old_major) > float(new_major)
            or float(old_minor) > float(new_minor)):
            print 'Error: Old Resolution bigger than new one!'

        fwhm_major = (math.sqrt(float(new_major) ** 2
                    - float(old_major) ** 2)
                    * const.arcsecInRad)
        fwhm_minor = (math.sqrt(float(new_minor) ** 2
                     - float(old_minor) ** 2)
                     * const.arcsecInRad)
        __gauss_smoothInit = open('gauss_smooth.init', 'w')
        __smooth_init_text = ('TASK\FILE "Input file" Y_NAME$ ' +
                          str(self.mapName) + '\n'
                         'TASK\FILE "Output smoothed image" X_NAME$ ' +
                         str(self.returnName(prefix=prefix,
                             resolution=new_resolution)) + '\n'
                         'TASK\REAL "Major axis of convolving gaussian" '
                         'MAJOR$ ' + str(fwhm_major) + '\n'
                         'TASK\REAL "Minor axis of convolving gaussian" '
                         'MINOR$ ' + str(fwhm_minor) + '\n'
                         'TASK\REAL "Position angle" PA$ ' + str(pa) + '\n'
                         'TASK\GO')

        __gauss_smoothInit.write(__smooth_init_text)
        __gauss_smoothInit.close()

        __convFile = open('temp.greg', 'w')
        __convFile.write('run gauss_smooth gauss_smooth.init /nowindow\n'
                            'exit')
        __convFile.close()
        os.system('greg -nw @temp.greg')
        os.system('rm temp.greg')
        return GildasMap(self.returnName(resolution=new_resolution))

    def quick_preview(self, save=False, filename=None):
        r"""
        Plotting the map and optionally save the figure.

        Parameters
        ----------

        save: True or False
            Choose wether or nor to save the figure.
        filename: string
            The filename to for the saved plot. If None defaults to 
            ``'quick_preview.eps'``.
        """
        pygreg.comm('dev im w')
        pygreg.comm('lim /rg')
        pygreg.comm('set box match')
        pygreg.comm('pl')
        pygreg.comm('box /abs')
        if save:
            filename = filename or 'quick_preview.eps'
            pygreg.comm('ha '+filename+' dev eps color')

    def toFits(self):
        r"""
        Converts the actual map to a Fits map.

        Returns
        -------

        FitsMap Object.

        Examples
        --------

        With:

        >>> map = gildasMap('M33_MIPS_24mum_JyB_5.gdf')
        >>> map = map.toFits()

        it is possible to continue working with the fits map, using the
        :class:`maps.fits.FitsMap` class.
        """
        comment = None
        for name in ['mean', 'velo', 'width']:
            if self.dataFormat == name:
                comment = [name]
        os.system('rm ' + self.returnName(dataFormat='fits', comments=comment))
        __convFile = open('temp.greg', 'w')
        __convFile.write('fits ' + self.returnName(dataFormat='fits',
                         comments=comment) + ' from ' + str(self.mapName) +
                         '\nexit\n')
        __convFile.close()
        os.system('greg -nw @temp.greg')
        self.fitsName = self.returnName(dataFormat='fits', comments=comment)
        return fits.FitsMap(self.fitsName)

    def toMiriad(self):
        r"""
        Converts the actual map to a Miriad map.

        Returns
        -------

        MiriadMap Object.

        Examples
        --------

        With:

        >>> map = gildasMap('M33_MIPS_24mum_JyB_5.gdf')
        >>> map = map.toMiriad()

        it is possible to continue working with the Miriad map, using
        :class:`maps.miriad.MiriadMap` class.
        """
        self.toFits()
        os.system('rm -rf ' + self.returnName(dataFormat=''))
        os.system('fits in=' + str(self.fitsName) + ' out=' +
                  self.returnName(dataFormat='') + ' op=xyin')
        self.miriadName = self.returnName(dataFormat='')
        self.dataFormat = ''
        return MiriadMap(self.miriadName)
