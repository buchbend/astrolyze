# Copyright (C) 2012, Christof Buchbender
# BSD Licencse
import math
import os
import string
import sys

from pysqlite2 import dbapi2 as sqlite
from scipy.ndimage import gaussian_filter
import pygreg

import main
import fits
import miriad

import astrolyze
from astrolyze.spectra import *
import astrolyze.functions.constants as const


class GildasMap(main.Map):
    r"""
    Wrapping GILDAS functions to use them inline with Python.
    """
    def __init__(self, map_name):
        r"""Initializes a Gildas map"""
        main.Map.__init__(self, map_name)
        self._init_map_to_greg()
        if self.dataFormat not in self.gildas_formats:
            print ('Exiting: Not a Gildas format (AFAIK). Supported'
                   'Formats Can be extended.')
            sys.exit()

    def _init_map_to_greg(self):
        pygreg.comm('image "' + self.map_name + '"')
        self._load_gildas_variables()
        pygreg.comm('let name "' +
                    self.map_name.replace('.' + self.dataFormat, '"'))
        pygreg.comm('let type "' + self.dataFormat + '"')

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
        self.ra_coordinate = self.vars.g_ra.__sicdata__
        self.dec_coordinate = self.vars.g_dec.__sicdata__
        # conversion to degrees
        self.ra_coordinate = self.ra_coordinate * const.r2d
        self.dec_coordinate = self.dec_coordinate * const.r2d
        self.central_coordinate_degrees = [self.ra_coordinate,
                                           self.dec_coordinate]
        self.central_coordinate_equatorial = astFunc.degrees_to_equatorial(
                                             self.central_coordinate_degrees)

    def set_defaults(self):
        r"""
        Reset all selection criteria.
        """
        pyclass.comm('set def')
        pyclass.comm('clear')

    def _regrid_to_arcsec(self, value):
        r"""
        Regrids the pixel size of the map to a multiple of arcseconds.

        Parameters
        ----------
        value : float
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
        self.cdelt1arcs = float(self.cdelt_1) / (value * (1. / 60 / 60))
        self.naxis1New = self.naxis_1 * math.sqrt(self.cdelt1arcs ** 2)
        self.crpix1New = self.naxis1New / (self.naxis_1 / self.crpix_1)
        self.comments += ['IP1']
        interpolate_init = open('interpolate.init', 'w')
        init_string = ('TASK\FILE \"Input file\" Y_NAME$ ' +
                         str(self.map_name) + '\n'
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
        interpolate_init.write(init_string)
        fileout.close()
        fileout = open('interpolate.greg', 'w')
        fileout.write('run interpolate interpolate.init /nowindow\nexit\n')
        fileout.close()
        os.system('greg -nw @interpolate.greg')

    def spectrum(self, coordinate, fileout=None, prefix=None,
                 create_spectrum=True):
        r""" Wrapper to the GILDAS/GREG spectrum command.

        Extracting a spectrum from a cube at a given positions. By default it
        also creates a 30m file readable by class from the table.

        Parameters
        ----------
        coordinate : list
            A list with the coordinates in floats in units of Degrees,
            or in string for equatorial coordinate.
        fileout : string
            The name of the table where the spectrum will be stored.
            Default is the same name as the map with ``".tab"`` as ending.
        prefix :
            The path to the folder where the newly created file will be
            stored.  Defaults to the prefix currently stored in self.prefix.
        create_spectrum : True or False
            Turn the creation of a 30m with the spectrum of ``False`` or on
            ``True``.

        Examples
        --------

        >>> from astrolyze.maps import *
        >>> map =  GildasMap('M33_PdBI_12co10_Tmb_22.0_2kms.gdf')
        >>> coordinate = ['1:34:7.00', '+30:47:52.00']
        >>> map.spectrum(coordinate)

        Creates M33_PdBI_12co10_Tmb_22.0_2kms.tab in the present folder.

        ..        Tested and working.
        """
        os.system('rm spectrum.init')
        if prefix is None:
            prefix = self.prefix
        if fileout is None:
            fileout = self.returnName(prefix=prefix, dataFormat='tab')
        else:
            fileout = prefix + fileout
        maskInit = open('spectrum.init', 'w')
        init_string = ('! Spectrum.init\n'
                       'TASK\FILE "Input file name" IN$ "' + self.map_name +
                       '"\nTASK\FILE "Output table name" OUT$ "'
                       + fileout
                       + '"\n'
                       'TASK\CHARACTER "Absolute coordinate" COORD$ "'
                       + coordinate[0] + ' ' + coordinate[1] + '"\n'
                       'TASK\GO')
        maskInit.write(init_string)
        maskInit.close()
        convFile = open('temp.greg', 'w')
        convFile.write('run spectrum spectrum.init /nowindow\nexit\n')
        convFile.close()
        os.system('greg -nw @temp.greg')
        os.system('rm temp.greg')
        os.system('rm spectrum.init')
        if create_spectrum:
            file_30m = self.returnName(comments=['extract'], prefix=prefix,
                                       dataFormat='30m')
            print 'file out ' + file_30m + ' single /over'
            print 'col x 1 y 2  /table ' + prefix + fileout
            pyclass.comm('file out ' + file_30m + ' single /over')
            pyclass.comm('col x 1 y 2  /table ' + prefix + fileout)
            pyclass.comm('model y x')
            pyclass.comm('write')

    def lmv(self, fileout=None, prefix=None):
        r"""Wrapper for the lmv command of GILDAS/CLASS.
        Extracts spectra from a spectral cube.

        Parameters
        ----------
        fileout : string
            The name of the class file to write the spectra to. Defaults to
            the map_name with .30m ending.

        prefix : string
            The path were the class file will be stores. Defaults to
            the current path.

        """
        if prefix is None:
            prefix = self.prefix
        if fileout is None:
            fileout = self.returnName(prefix=prefix, dataFormat='30m')
        else:
            fileout = prefix + fileout
        # pyclass.comm('file out {} single /over'.format(fileout))
        convFile = open('temp.greg', 'w')
        convFile.write('file out {} single /over\n'.format(fileout))
        convFile.write('lmv {}\nexit\n'.format(self.map_name))
        convFile.close()
        os.system('class -nw @temp.greg')
        os.system('rm temp.greg')
        # pyclass.comm('lmv {}'.format(self.map_name))
        return ClassSpectra(fileout)

    def mask(self, polygon, prefix=None):
        r"""
        Wrapper for the GREG task mask:

        Parameters
        ----------

        polygon : string
            path to a GILDAS polygon file with ending ``".pol"``

        prefix : string
            The path where the output is to be stored if different
            from the current prefix of the map.

        Returns
        -------

        mapObject : The masked map object.

        Examples
        --------

        >>> map.mask('poly/sourceA.pol')
        """
        if prefix is None:
            prefix = self.prefix
        maskInit = open('mask.init', 'w')
        comment = [polygon.split('/')[-1].replace('.pol', '')]
        init_string = ('TASK\FILE "Polygon data file" POLYGON$ "' +
                         str(polygon) + '" \n'
                         'TASK\FILE "Input image" Y_NAME$ "' +
                         str(self.map_name) + '"\n'
                         'TASK\FILE "Output image" X_NAME$ "'
                         + str(self.returnName(prefix=prefix,
                                               comments=comment)) + '"\n'
                         'TASK\LOGICAL "Mask inside (.true.) or outside '
                         '(.false.)" MASK_IN$ NO\n'
                         'TASK\LOGICAL "Modify blanking value" MODIFY$ NO\n'
                         'SAY "If YES, Fill information Below"\n'
                         'TASK\REAL "New blanking value" BLANKING$ 0 \n'
                         'TASK\GO')
        maskInit.write(init_string)
        maskInit.close()
        os.system('more mask.init')
        convFile = open('temp.greg', 'w')
        convFile.write('run mask mask.init /nowindow\nexit\n')
        convFile.close()
        os.system('greg -nw @temp.greg')
        os.system('rm temp.greg')
        os.system('rm mask.init')
        return GildasMap(self.returnName(prefix=prefix,  comments=comment))

    def reproject(self, template=None, coord=None, prefix=None,
                  keep_pixsize=False):
        r"""
        Wraps the GREG task reproject. Either use template *or* coord.

        Parameters
        ----------
        template : string
            Full path to a map in GDF Format whose central
            coordinate and pixel size will serve as a template.
        coord : list
            List of coordinate strings in RA DEC (J2000) that
            will become the new centre of the map.
        prefix : string
            The path where the output is to be stored if different
            from the current prefix of the map. If None the current
            self.prefix of the GildasMap instance is used.
        keep_pixsize : bool
            If False reproject guesses the new pixel sizes after reprojection
            these are normally smaller than the original ones.
            If True the old pixel sizes are enforced.

        Returns
        -------
        GildasMap Object : Instance for the reprojected map.

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
        reproInit = open('reproject.init', 'w')
        if prefix is None:
            prefix = self.prefix
        if template:
            os.system('cp {} template.gdf'.format(template))
            template = 'template.gdf'
            comment = ['repro']
            init_string = ('TASK\FILE "Input file" Z_NAME$ "' +
                             str(self.map_name) + '"\n'
                             'TASK\FILE "Output file" X_NAME$ "' +
                             str(self.returnName(prefix=prefix,
                                                 comments=comment))
                             + '"\n'
                             'TASK\LOGICAL "Reproject on a third file\'s '
                             'projection [YES|NO]" TEMPLATE$ YES\n'
                             'TASK\FILE "This Template File Name" Y_NAME$ "' +
                             str(template) + '"\n'
                             'TASK\CHARACTER "Projection type" PROJECTION$ '
                             '"RADIO" /CHOICE     NONE GNOMONIC ORTHOGRAPHIC  '
                             'AZIMUTHAL STEREOGRAPHIC LAMBERT AITOFF RADIO\n'
                             'TASK\CHARACTER "Coord.System '
                             '[EQUATORIAL [epoch]|'
                             'GALACTIC]" SYSTEM$ "EQUATORIAL 2000"\n'
                             'TASK\CHARACTER "Coord. 1 of projection center. '
                             'Could be UNCHANGED"     CENTER_1$ "01:34:11.8" '
                             '/CHOICE UNCHANGED *\n'
                             'TASK\CHARACTER "Coord. 2 of projection '
                             'center" CENTER_2$ "+30:50:23.4" '
                             '/CHOICE UNCHANGED '
                             '*\n'
                             'TASK\REAL "Position angle of projection" ANGLE$ '
                             '0.000000000000000 /RANGE 0 360\n'
                             'TASK\INTEGER "Dimensions of output image '
                             '[0 0 mean guess]" DIMENSIONS$[2]  0 0\n'
                             'TASK\REAL "First axis conversion formula [0 0 0 '
                             'mean guess]" AXIS_1$[3]  0 0 0\nTASK'
                             '\REAL "Second '
                             'axis conversion formula [0 0 0 mean guess]" '
                             'AXIS_2$[3] 0 0 0\n'
                             'TASK\LOGICAL "Change blanking value" '
                             'CHANGE$ NO\n'
                             'TASK\REAL "New blanking value and tolerance" '
                             'BLANKING$[2]  0 0\nTASK\GO\n')
            reproInit.write(init_string)
        if coord:
            comment = ['repro']
            init_string = ('TASK\FILE "Input file" Z_NAME$ "' +
                             str(self.map_name) + '"\n'
                             'TASK\FILE "Output file" X_NAME$ "' +
                             str(self.returnName(prefix=prefix,
                                                 comments=comment))
                             + '"\n'
                             'TASK\LOGICAL "Reproject on a third file\'s '
                             'projection [YES|NO]" TEMPLATE$ NO\n'
                             'TASK\FILE "This Template File Name" Y_NAME$ ""\n'
                             'TASK\CHARACTER "Projection type" PROJECTION$ '
                             '"RADIO" /CHOICE     NONE GNOMONIC ORTHOGRAPHIC  '
                             'AZIMUTHAL STEREOGRAPHIC LAMBERT AITOFF RADIO\n'
                             'TASK\CHARACTER "Coord.System '
                             '[EQUATORIAL [epoch]|'
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
                init_string += ('TASK\INTEGER "Dimensions of output image [0'
                                  '0 mean guess]" DIMENSIONS$[2] ' +
                                  str(self.naxis_1) + ' ' + str(self.naxis_1) +
                                  '\n' 'TASK\REAL "First axis conversion'
                                  'formula [0 0 0 ' 'mean guess]" AXIS_1$[3] '
                                  + str(self.crpix_1) + ' '
                                  + str(self.crval_1) + ' '
                                  + str(self.cdelt_1) + '\n' 'TASK\REAL'
                                  '"Second axis conversion formula [0 0 0 '
                                  'mean guess]" AXIS_2$[3] ' +
                                  str(self.crpix_2) + ' ' +
                                  str(self.crval_2) + ' ' +
                                  str(self.cdelt_2) + '\n')
            if not keep_pixsize:
                init_string += ('TASK\INTEGER "Dimensions of output image [0'
                                  '0 mean guess]" DIMENSIONS$[2] 0 0 \n'
                                  'TASK\REAL "First axis conversion formula '
                                  '[0 0 0 mean guess]" AXIS_1$[3] 0 0 0 \n '
                                  'TASK\REAL "Second axis conversion formula '
                                  '[0 0 0 mean guess]" AXIS_2$[3] 0 0 0 \n')
            init_string += ('TASK\LOGICAL "Change blanking value" CHANGE$ '
                              'YES\n TASK\REAL "New blanking value and '
                              'tolerance" BLANKING$[2]  0 0\n'
                              'TASK\GO\n')
            reproInit.write(init_string)
        reproInit.close()
        # write the Greg Script
        convFile = open('temp.greg', 'w')
        convFile.write('run reproject reproject.init /nowindow\nexit\n')
        convFile.close()
        os.system('more reproject.init')
        # execute the Greg Script and then delete it
        os.system('greg -nw @temp.greg')
        os.system('rm temp.greg')
        os.system('rm reproject.init')
        os.system('rm template.gdf')
        return GildasMap(self.returnName(prefix=prefix, comments=comment))

    def moments(self, velo_range=[0, 0], threshold=0,
                smooth='YES', prefix=None, comment=None):
        r""" Wraps the GILDAS/GREG task moments.

        Creates the first three moments of the map.

        Parameters
        ----------
        velo_range : list
            Velocity range for the integration.
        threshold : float
            Value under which pixels are blanked.
        smooth : string
            One of Either ``"NO"`` or ``"YES"``. Controls
            if the map is smoothed in velocity before applying the cut
            threshold. Getting rid of noise peaks over the threshold.
            Defaults to ``'YES'``
        prefix : string
            The path where the output is to be stored if different
            from the current prefix of the map.
        comment : string
            Optional comments to be added to the new map name.

        Returns
        -------
        mean : MapObject
            The zeroth moment, i.e. the integrated intensity, is returned as a
            GildasMap object.
        """
        print self.map_name
        if comment is None:
            comment = []
        save_comments = self.comments
        self.comments = []
        if prefix is None:
            prefix = self.prefix
        momentInit = open('moments.init', 'w')
        init_string = ('TASK\FILE "Input file name" IN$ "' +
                         str(self.map_name) + '"\n'
                         'TASK\FILE "Output files name (no extension)" '
                         'OUT$ "' +
                         self.returnName(resolution='dummy', prefix=prefix,
                                         comments=[]).replace('.' +
                                                              self.dataFormat,
                                                              '') + '"\n'
                         'TASK\REAL "Velocity range" VELOCITY$[2]  ' +
                         str(velo_range[0]) + ' ' + str(velo_range[1]) + '\n'
                         'TASK\REAL "Detection threshold" THRESHOLD$ ' +
                         str(threshold) + ' \n'
                         'TASK\LOGICAL "Smooth before detetction ?" SMOOTH$ ' +
                         str(smooth) + ' \n'
                         'TASK\GO\n')
        momentInit.write(init_string)
        momentInit.close()
        convFile = open('temp.greg', 'w')
        convFile.write('run moments moments.init /nowindow\nexit\n')
        convFile.close()
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
        os.system('rm moments.init')
        print self.returnName(prefix=prefix,
                              comments=save_comments,
                              dataFormat='mean')
        return GildasMap(self.returnName(prefix=prefix,
                                         comments=save_comments,
                                         dataFormat='mean'))

    def goRot(self, angle, prefix=None):
        r"""
        Wrapper for the GREG go rot command, which rotates maps around their
        central coordinate stored in the header.

        Parameters
        ----------
        angle : float [deg]
            Rotation angle.
        prefix : string
            The path where the output is to be stored if different
            from the current prefix of the map.

        Returns
        -------
        GildasMap Object : Instance for the reprojected map.

        Examples
        --------
        >>> map.goRot(45)

        To change the central coordinate first use
        :py:func:`maps.gildas.GildasMap.reproject` e.g.:

        >>> map = map.reproject(coord=['new_RA_string','new_DEC_string'])
        >>> map.goRot(45)

        """
        rotate = open('rotateTemp.greg', 'w')
        rotate.write('let name ' +
                       self.map_name.replace('.' + self.dataFormat, '') +
                       '\n let type ' + str(self.dataFormat) + '\n'
                       'let angle ' + str(angle) + '\n'
                       'go rot\n'
                       'exit\n')
        rotate.close()
        os.system('greg -nw @rotateTemp.greg')
        if len(str(angle).split('.')) > 1:
            comments += [('rot' + str(str(angle).split('.')[0]) +
                          'p' + str(str(angle).split('.')[1]) + 'deg')]
        else:
            comments += ['rot' + str(angle) + 'deg']

        os.system('mv ' + self.map_name.replace('.'+ self.dataFormat, '' ) +
                  '-rot' + str(angle) + 'deg.'
                  + str(self.dataFormat) + ' '
                  + str(self.returnName(prefix=prefix, comments=comment)))
        os.system('rm -rf rotateTemp.greg')
        return GildasMap(self.returnName(prefix=prefix, comments=comment))

    def slice(self, coordinate1, coordinate2, prefix=None, comment=None):
        r""" Wrapper for the GREG task slice. Producing Position velocity cuts
        trough a map between coordinate1 and coordinate2.

        Parameters
        ----------

        coordinate1 : string
            The coordinate where the cut trough the map starts.
        coordinate2 : string
            The coordinate where the cut trough the map ends.

        Returns
        -------

        A GildasMap object containing the slice.

        Notes
        -----

        This only works with cubes.
        """
        if prefix is None:
            prefix = self.prefix
        if comment is None:
            comment = ['sliced']
        elif comment is not None and '[' not in comment:
            comment = [comment]
        sliceInit = open('slice.init', 'w')
        init_string = ('!\n! slice.init\n' +
                         'TASK\FILE "Input file name" INPUT_MAP$ "' +
                         str(self.map_name) + '"\n' +
                         'TASK\FILE "Output file name" OUTPUT_MAP$ "' +
                         str(self.returnName(prefix=prefix,
                                             comments=comment)) + '"\n' +
                         'TASK\CHAR "Starting point" START$ "' +
                         str(coordinate1) + '"\n' +
                         'TASK\CHAR "Ending point" END$ "' +
                         str(coordinate2) + '"\n' +
                         'TASK\GO')
        sliceInit.write(init_string)
        sliceInit.close()
        convFile = open('temp.greg', 'w')
        convFile.write('run slice slice.init /nowindow\nexit\n')
        convFile.close()
        os.system('greg -nw @temp.greg')
        os.system('rm temp.greg')
        os.system('rm slice.init')
        return GildasMap(self.returnName(prefix=prefix, comments=comment))

    def smooth(self, new_resolution, old_resolution=None, prefix=None):
        r"""Wrapper for the GILDAS/GREG task gauss_smooth.

        Parameters
        ----------
        new_resolution : float or list
            The resulting resolution after the smoothing.
            It can be:

            1. a float: i.e. the final major and minor beamsize.
               The position angle will default to 0.
            2. a list with two floats: [major_axis, minor_axis]. The
               position angle defaults to 0.
            3. a list with three floats: [major_axis, minor_axis,
               position_angle].

        old_resolution : float or list
            Same format as new_resolution. Defaults to self.resolution of the
            map instance.
        prefix : string
            The path where the output is to be stored if different
            from the current prefix of the map.

        Notes
        -----

        .. warning::

            The gauss_smooth Task from GILDAS only gives correct output units
            when the map is on a temperature or \"per pixel\" scale.  **Maps in
            Jy/Beam won't be in Jy/Beam after smoothing.**

        """
        if prefix is None:
            prefix = self.prefix
        if old_resolution is None:
            old_major = self.resolution[0]
            old_minor = self.resolution[1]
            pa = self.resolution[2]
        if old_resolution is not None:
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

        if ((float(old_major) > float(new_major)
             or float(old_minor) > float(new_minor))):
            self.log.error('Old Resolution bigger than new one!')

        fwhm_major = (math.sqrt(float(new_major) ** 2
                                - float(old_major) ** 2)
                      * const.arcsecInRad)
        fwhm_minor = (math.sqrt(float(new_minor) ** 2
                                - float(old_minor) ** 2)
                      * const.arcsecInRad)
        smooth_init_text = ('TASK\FILE "Input file" Y_NAME$ ' +
                              str(self.map_name) + '\n'
                              'TASK\FILE "Output smoothed image" X_NAME$ ' +
                              str(self.returnName(prefix=prefix,
                                                  resolution=new_resolution)) +
                              '\n TASK\REAL "Major axis of '
                              'convolving gaussian" '
                              'MAJOR$ ' + str(fwhm_major) + '\n'
                              'TASK\REAL "Minor axis of convolving gaussian" '
                              'MINOR$ ' + str(fwhm_minor) + '\n'
                              'TASK\REAL "Position angle" PA$ ' + str(pa) +
                              '\n TASK\GO')
        with open('gauss_smooth.init', 'w') as output:
            output.write(smooth_init_text)

        convFile = open('temp.greg', 'w')
        convFile.write('run gauss_smooth gauss_smooth.init /nowindow\n'
                         'exit')
        convFile.close()
        os.system('greg -nw @temp.greg')
        os.system('rm temp.greg')
        os.system('rm gauss_smooth.init')
        return GildasMap(self.returnName(resolution=new_resolution))

    def custom_go_spectrum(self, coordinate=False, size=False, angle=0):
        r"""
        This function uses the ``go spectrum`` command from GREG to plot
        the spectra in a region given by ``size`` around the cooridinate given
        by ``coordinate``.

        Parameters
        ----------
        coordinate : list
            A list with the coordinates in floats in units of Degrees, or in
            string for equatorial coordinates. Default to ``False`` which means
            that the center of the map, determined from the map header, is
            used.
        size : list
            The region around the ``coordinate`` from which spectra are plotted
            in arcsec, e.g. size = [50, 50] means a region of 50x50 arcsec
            around the given cooridnate. Defaults to None, which translates to
            size = [0, 0] which in turn is interpreted as  the full map size
            by GREG.
        angle : float [degrees]
            Needed if the map is rotated to get the correct offsets. Defaults
            to 0.
        """
        # go spectrum reads the variable center which gives the coordinates
        # in offsets from the central coordinate in the header
        if coordinate:
            offset = astFunc.calc_offset(self.central_coordinate_equatorial,
                                         coordinate, angle,
                                         output_unit='arcsec')
        else:
            offset = [0, 0]
        print offset
        if not size:
            size = [0, 0]
        comment = ['go-spectrum']
        name = self.returnName(comments=comment, dataFormat='eps')
        string = ('let center ' + str(offset[0]) + ' ' + str(offset[1]) + '\n'
                  'let size ' + str(size[0]) + ' ' + str(size[1]) + '\n'
                  'let name ' + self.map_name.replace('.' +
                                                      self.dataFormat,
                                                      '') + '\n'
                  'let type ' + self.dataFormat + '\n'
                  'go spectrum \n'
                  'ha ' + name + ' /dev eps color /over \n'
                  'exit \n'
                  )
        print string
        fileout = open('temp.greg', 'w')
        fileout.write(string)
        fileout.close()
        os.system('greg -nw @temp.greg')
        os.system('rm -rf temp.greg')

    def save_figure(self, name=None):
        r"""
        Helper function that saves the current plot.
        """
        name = name or self.returnName(dataFormat='eps')
        pygreg.comm('ha ' + name + '/dev eps color')

    def quick_preview(self, save=False, filename=None, window=True,):
        r"""
        Plotting the map and optionally save the figure.

        Parameters
        ----------
        save : True or False
            Choose wether or nor to save the figure.
        filename : string
            The filename to for the saved plot. If None defaults to
            ``'quick_preview.eps'``.
        window : True or False
            Choose whether the image display is opened or not.
            Default True.
        """
        pygreg.comm('clear')
        pygreg.comm('image ' + self.map_name)
        if window:
            pygreg.comm('dev im w')
        pygreg.comm('lim /rg')
        pygreg.comm('set box match')
        pygreg.comm('greg2\\pl')
        pygreg.comm('greg\\box /abs')
        if save:
            filename = filename or 'quick_preview.eps'
            if '.eps' in filename:
                pygreg.comm('ha ' + filename + ' /dev eps color /over')
            if '.png' in filename:
                pygreg.comm('ha ' + filename + ' /dev png color /over')


    def get_spectra_from_cube(self, coordinate, angle=0, prefix=None,
                              accuracy=2):
        spectra = self.lmv()
        spectra = spectra.get_spectra_from_cube(coordinate, angle=0,
                                                prefix=None, accuracy=2)
        return spectra

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
        convFile = open('temp.greg', 'w')
        convFile.write('fits ' + self.returnName(dataFormat='fits',
                         comments=comment) + ' from ' + str(self.map_name) +
                         '\nexit\n')
        convFile.close()
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
