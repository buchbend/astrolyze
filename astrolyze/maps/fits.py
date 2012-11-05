# Copyright (C) 2009-2012, Christof Buchbender
# BSD License (License.txt)
import os
import string
import sys
import pyfits
import pywcs
import math

from pysqlite2 import dbapi2 as sqlite
from scipy.ndimage import gaussian_filter

import numpy as np

import main
import gildas
import miriad

import astrolyze.functions.constants as const
from astrolyze.functions import astro_functions as astFunc
from astrolyze.functions import units


class FitsMap(main.Map):
    r"""
    Fits Map manipulation making extensive use of the
    pyfits package.
    """
    def __init__(self, map_name, name_convention=True):
        main.Map.__init__(self, map_name, name_convention)
        if self.dataFormat not in self.fits_formats:
            print 'Exiting: not the right format'
            sys.exit()
        self.__get_data()
        self.__get_header()
        self.wcs = pywcs.WCS(self.hdulist[0].header)
        self.headerKeywords = {}
        self.headerKeywords['source'] = ['OBJECT', 'SOURCE', 'TARNAME']
        self.headerKeywords['telescope'] = ['TELESCOP']
        self.headerKeywords['species'] = ['LINE', 'FREQUENCY']
        self.headerKeywords['fluxUnit'] = ['BUNIT']
        self.headerKeywords['resolution'] = ['BMAJ', 'BMIN']
        self.get_pixel_size()

    def __get_data(self):
        r"""
        Creates a self.data variable containing an array with the maps data.
        This function is automatically executed every time a Fits map is opened
        """
        try:
            self.data = self.hdulist[0].data
        except:
            self.hdulist = pyfits.open(self.map_name)
            self.data = self.hdulist[0].data

    def __get_header(self):
        r"""
        Creates a self.header variable containing a dictionary with the header
        information of the map. This function is automatically executed every
        time a Fits map is opened.
        """
        try:
            self.header = self.hdulist[0].header
        except:
            self.hdulist = pyfits.open(self.map_name)
            self.header = self.hdulist[0].header

    def update_file(self, backup=False):
        r"""
        Writing changes to the self.data and/or self.header to the current
        file.

        Parameters
        ----------
        backup : True or False
            If True a copy of the original file is created having the extension
            ``"_old"`` after the file endind, i.e. some_name.fits -> some_name.
            fits_old.

        Returns
        --------
        FitsMap : Instance

        Notes
        -----
        If all variables that define the map name () are unchanged the current
        file is overwritten else

        Examples
        --------
        """
        map_name = self.returnName()
        if backup is True:
            os.system('cp ' + str(map_name) + ' ' + str(map_name) + '_old')
        os.system('rm ' + str(map_name))
        try:
            pyfits.writeto(map_name, self.data, self.header)
        except:
            print 'Problem with the header or data format. -> Exit!'
            raise SystemExit
        return FitsMap(map_name)

    def updateHeader(self):
        try:
            self.header
        except:
            self.get_header()
        self.changemap_name()
        # Print source
        updateCheck = 0
        for i in self.headerKeyWords['source']:
            if str(i) in self.header.keys():
                self.header.update(str(i), self.source)
                updateCheck = 1
        if updateCheck == 0:
            self.header.update(str(self.headerKeyWords['source'][0]),
                               str(self.source))
        # Print telescope
        updateCheck = 0
        for i in self.headerKeyWords['telescope']:
            if str(i) in self.header.keys():
                self.header.update(str(i), self.telescope)
                updateCheck = 1
        if updateCheck == 0:
            self.header.update(str(self.headerKeyWords['telescope'][0]),
                               str(self.telescope))
        # Print species
        updateCheck = 0
        for i in self.headerKeyWords['species']:
            if str(i) in self.header.keys():
                self.header.update(str(i), self.species)
                updateCheck = 1
        if updateCheck == 0:
            self.header.update(str(self.headerKeyWords['species'][0]),
                               str(self.species))
        # Update fluxUnit
        for i in self.headerKeyWords['fluxUnit']:
            if str(i) in self.header.keys():
                self.header.update(str(i), self.fluxUnit)
                updateCheck = 1
        if updateCheck == 0:
            self.header.update(str(self.headerKeyWords['fluxUnit'][0]),
                               str(self.fluxUnit))
        # Update resolution
        for i in self.headerKeyWords['resolution']:
            if str(i) in self.header.keys():
                self.header.update(str(i),
                                   str(float(self.resolution[0]) / 60 / 60))
                updateCheck = 1
        if updateCheck == 0:
            self.header.update(str(self.headerKeyWords['resolution'][0]),
                               str(self.resolution[0]))
        # update header information about the minimum and maximum of the map
        try:
            self.header['DATAMIN'] = self.data.min() - 0.1 * self.data.min()
        except:
            self.header.update('DATAMIN', self.data.min() - 0.1 *
                               self.data.min())
        try:
            self.header['DATAMAX'] = self.data.max() + 0.1 * self.data.max()
        except:
            self.header.update('DATAMAX', self.data.max() + 0.1 *
                               self.data.max())
        os.system('cp ' + str(self.map_name) + ' ' +
                  str(self.map_name) + '_old')
        os.system('rm ' + str(self.map_name))
        pyfits.writeto(self.map_name, self.data, self.header)

    def get_pixel_size(self):
        r""" Calculates the Area of a pixel in m^2 and steradians if distance
        is given. If not only  steradians are calculated.
        """
        _pixSize = np.asarray([float(math.fabs(self.header['CDELT1']) *
                                     const.d2r),
                               float(math.fabs(self.header['CDELT2'])
                                     * const.d2r)])
        _pixSizeSterad = _pixSize[0] * _pixSize[1]
        self.pixelSizeSterad = math.sqrt((float(_pixSizeSterad)) ** 2)
        if self.distance is not None:
            _pixSize = _pixSize * self.distance * const.pcInM
            _pixSize = _pixSize[0] * _pixSize[1]
            self.pixelSizeM2 = math.sqrt((float(_pixSize)) ** 2)

    def __smooth(self, newRes, oldRes=None, scale='0.0'):
        # TODO: In Development. It's not sure if it ever will be finished.
        # No need to reinvent the wheel!
        r"""
        ..  warning::

            Still in Development. Maybe Never finished.
            **Do Not Use**
        """
        if oldRes is None:
            oldMajor = self.resolution[0]
            oldMinor = self.resolution[1]
            pa = self.resolution[2]
        if float(oldRes) > float(newRes):
            print 'Error old Resolution bigger than new one!'
        # calculate the fwhm for the convolving gaussian
        newRes = newRes * 4.848e-6
        oldRes = oldRes * 4.848e-6
        fwhmMajor = math.sqrt(float(newRes) ** 2 - float(oldMajor) ** 2)
        fwhmMinor = math.sqrt(float(newRes) ** 2 - float(oldMinor) ** 2)
        os.system('rm -rf ' + str(self.returnName(resolution=str(newRes))))
        # Have to check for Scaling!!!!!!
        while len(self.data) == 1:
            self.data = self.data[0]
        self.data = gaussian_filter(self.data, fwhm)
        # set the map_name and the resolution to the new values<
        self.map_name = self.returnName(resolution=str(newRes))
        self.resolution = newRes

    def _cut_map(self, x1y1, x2y2, pix_or_coord='coord'):
        # TODO: Check if still valid!!!!
        r"""
        Cutting an rectangle out of a map. Giving the corners in
        coordinates or in pixels.

        Parameters
        ----------
        x1y1 : list
            The upper right corner of the rectangle to cut out.
            Either in:
                * pixel coordinated [x1, y1]
            Or:
                * equatorial coordinates ['RA','DEC']
        x2y2 : list
            The lower left corner in the same format as x1y1.
        pix_or_coord : string
            Either ``"pix"`` or ``"coord"`` choosing what **x1y1**
            and **x2y2** represents.

        Notes
        -----
        This procedure cuts only rectangles paralell to the sides of the map.

        ..  warning::
            Old function. Functionality not guaranteed. maybe not really
            useful. Test or remove.
        """
        if  pix_or_coord == 'coord':
            x1y1 = self.sky2pix(x1y1)
            print 'x1y1', x1y1
            x2y2 = self.sky2pix(x2y2)
            print 'x2y2', x2y2
            if (x1y1 or x2y2) == 'Error':
                return 'Coordinates are not correct'
        if pix_or_coord == 'pix':
            pass
        print x1y1, x2y2
        x1 = math.floor(float(x1y1[0]))
        x2 = math.floor(float(x2y2[0]))
        y1 = math.floor(float(x1y1[1]))
        y2 = math.floor(float(x2y2[1]))
        try:
            self.data
        except:
            self.get_data()
        # needed if to check the dimension of the data array
        if len(self.data.shape) > 2:
            for i in range(len(self.data.shape) - 2):
                self.data = self.data[0]
        print x1, x2, y1, y2
        self.cutData = self.data[x1:x2, y1:y2]
        self.NewCRVal1 = ((x1 - math.floor(self.header['CRPIX1']))
                          * self.header['CDELT1'] + self.header['CRVAL1'])
        self.header['CRVAL1'] = self.NewCRVal1
        self.NewCRVal2 = ((y1 - math.floor(self.header['CRPIX2']))
                          * self.header['CDELT2'] + self.header['CRVAL2'])
        self.header['CRVAL2'] = self.NewCRVal2
        print (str(len(range(int(x1), int(x2)))),
               str(len(range(int(y1), int(y2)))))
        self.header['NAXIS1'] = len(range(int(x1), int(x2)))
        self.header['NAXIS2'] = len(range(int(y1), int(y2)))
        self.header['CRPIX2'] = 1
        self.header['CRPIX2'] = 1
        print self.header
        self.data = self.cutData
        self.comments += ['cut']
        self.update_file()

    def _strip(self, coords, radial=None, centerCoord=None, distance=None):
        # TODO: OLD code. Check or remove.
        r"""
        Extracts a linear cut trough a map between two coordinates.

        Distance in kpc.

        Notes
        -----
        ..  warning:

            Old code. Test and document or remove.

        """
        while len(self.data[0]) == 1:
            self.data = self.data[0]
        x1, y1 = self.sky2pix(coords[0])
        x2, y2 = self.sky2pix(coords[1])
        if centerCoord is not None:
            centx, centy = self.sky2pix(centerCoord)
            print centy, centx
        values = [[], []]
        # in arcsec
        self.pixSize1 = float(self.header['CDELT1']) / (1. / 60 / 60)
        self.pixSize2 = float(self.header['CDELT2']) / (1. / 60 / 60)
        if x1 != x2 and y1 != y2:
            print ('This version of strip only supports linear cuts.'
                   'Please rotate your map.')
        if x1 == x2:
            if y1 > y2:
                loop = range(y2, y1)
            if y1 < y2:
                loop = range(y1, y2)
            for i in loop:
                values[0] += [self.data[i][x1]]
                if radial == 1:
                    values[1] += [(i - centy) * self.pixSize2 / 60]
                else:
                    values[1] += [i]
        if y1 == y2:
            if x1 > x2:
                loop = range(x2, x1)
            if x1 < x2:
                loop = range(x1, x2)
            for i in loop:
                values[0] += [self.data[y1][i]]
                if radial is True:
                    values[1] += [(i - centx) * self.pixSize1 / 60]
                else:
                    values[1] += [i]
        return values

    def read_flux(self, position):
        r"""
        Returns the value of the pixel that corresponds to the
        given positions of RA, DEC (J2000) in units of equatorial
        coordinates or degrees.

        Parameters
        ----------

        position : list
            The position in RA,DEC where the aperture is to be applied.
            The Format has to be either:

                * ['RA','DEC'] with strings representing equatorial
                  coordinates, e.g. ['01:34:32.8', '+30:47:00.6'].
            or:
                * [RA, DEC] where RA and DEC being the coordinates in Grad.

        Returns
        -------

        flux : float
            The flux at the given position.

        See Also
        --------

        sky2pix, astFunc.equatorial_to_degrees, wcs.wcs_sky2pix

        """
        x, y = self.sky2pix(position)
        return self.data[y][x]

    def read_aperture(self, position, apertureSize=0, backgroundSize=0,
                      output=False, annotation=False, newAnnotation=False):
        # TODO: needs better scheme for the header keywords. Migth not work
        # with all maps.
        r"""
        Extract the sum and mean flux inside an aperture of a given size
        and at a given position..

        This function can be used to read the flux in the area of a circular
        aperture, as well as to correct for the background flux.

        Parameters
        ----------
        position : list
            The position in RA,DEC where the aperture is to be applied.
            The Format has to be either:

                * ['RA','DEC'] with strings representing equatorial
                  coordinates, e.g. ['01:34:32.8', '+30:47:00.6'].
            or:
                * [RA, DEC] where RA and DEC being the coordinates in Grad.

        apertureSize : float [arcsec]
            The diameter of the aperture to be applied.

        backgroundSize : float [arcsec]
            The Size of the Anulli in which the background is to be
            estimated. The number to be given here correspond to the diameter
            of the circle in [arcsec ] descibing the outer border of the
            annuli, measured from the position given in position. Thus, the
            background is measurd in the ring described by apertureSize and
            backgroundSize. Default is 0 and thus **no background substaction**
            is applied.

        output : True or False
            If True the function reports the calculated values during
            execution.

        annotion : logical
            If True a kvis annotation file ``"apertures.ann"`` containing the
            aperture used to integrate the flux is created. Default is False,
            i.e. not to create the aperture.

        newAnnotation : logical
            If True ``"apertures.ann"`` is overwritten. If False an old
            ``"apertures.ann"`` is used to append the new apertures. If it not
            exists a new one is created. The latter is the default.

        Returns
        -------
        List :  [Sum, Mean, Number of pixels]

        Notes
        -----

        The pixel sizes have to be quadratic for the algorithm to work. It
        measures a circle by counting the pixels from the central pixel
        corresponding to the given coordinate.

        """
        # Calculate the x,y Pixel coordinate of the position.
        xCenter, yCenter = self.sky2pix(position)
        # Read the pixel Size (normaly in Grad FITS-Standard)
        try:
            self.pixSize1 = float(self.header['CDELT1'])
            self.pixSize2 = float(self.header['CDELT2'])
        except:
            # sometimes the header keywords are different.
            try:
                self.pixSize1 = float(self.header['XPIXSIZE'])
                self.pixSize2 = float(self.header['YPIXSIZE'])
                self.pixSize1 = self.pixSize1 / 60 / 60
                self.pixSize2 = self.pixSize2 / 60 / 60
            except:
                print 'No Header keyword for the pixsize found'
                sys.exit()
         # Check if the pixel size is quadratic, if not the function does not
         # work correctly (yet)
        if ((((self.pixSize1) ** 2 - (self.pixSize2) ** 2) ** 2 <
             ((self.pixSize1) ** 2) * 0.01)):
            if output is True:
                print ('Pixel Size X: ' + str(self.pixSize1 / (1. / 60 / 60))
                       + ' arcsec')
                print ('Pixel Size Y: ' + str(self.pixSize2 / (1. / 60 / 60))
                       + ' arcsec')
                print ('Pixel Size quadratic to 1%, Go on with '
                       'calculation of Apertures')
        else:
            if output is True:
                print self.pixSize1, self.pixSize2
                print ('Pixel Size *NOT* quadratic: Exit the program, '
                       'please regrid!')
            sys.exit()
        # Convert the apertureSizes to Grad f
        apertureSizeInGrad = (1. / 60 / 60) * (apertureSize) / 2
        # Calculate how many pixels are needed to cover the radius of the
        # aperture
        apertureSteps = math.sqrt(math.floor(apertureSizeInGrad /
                                             self.pixSize1) ** 2)
        # same for the backround annuli
        if backgroundSize != 0:
            backgroundSizeInGrad = (1. / 60 / 60) * (backgroundSize) / 2
            backgroundSteps = math.sqrt(math.floor(backgroundSizeInGrad /
                                                   self.pixSize1) ** 2)
            x_max, x_min = xCenter + backgroundSteps, xCenter - backgroundSteps
            y_max, y_min = yCenter + backgroundSteps, yCenter - backgroundSteps
        else:
            x_max, x_min = xCenter + apertureSteps, xCenter - apertureSteps
            y_max, y_min = yCenter + apertureSteps, yCenter - apertureSteps
        # calculate the x and y range in Pixel coordinates that are needed at
        # most for the calculation of the aperture and background
        # Now: Initialize variables
        sum = 0  # will store the sum of all Pixels in the aperture
        N = 0   # counts the number of pixels in the aperture to be able to
                # calculate the mean
        sumBackgnd = 0  # the sum of all pixels in the background
        NBackgnd = 0   # number of pixels in the background
        # calculate the background
        if backgroundSize != 0:
            for x in range(int(x_min), int(x_max)):
                for y in range(int(y_min), int(y_max)):
                    if ((math.floor(math.sqrt((xCenter - x) ** 2 +
                                              (yCenter - y) ** 2)) <
                         backgroundSteps and math.floor(math.sqrt((xCenter - x)
                                                                  ** 2 +
                                                                  (yCenter - y)
                                                                  ** 2))
                         > apertureSteps)):
                        sumBackgnd += self.data[y][x]
                        NBackgnd += 1
            backgndMean = sumBackgnd / NBackgnd
            print backgndMean
        # caclulate the sum in the aperture
        for x in range(int(x_min), int(x_max)):
            for y in range(int(y_min), int(y_max)):
                if (math.floor(math.sqrt((xCenter - x) ** 2 + (yCenter - y) **
                                         2)) < apertureSteps):
                    if backgroundSize != 0:
                        sum += self.data[y][x] - backgndMean
                    else:
                        sum += self.data[y][x]
                    N += 1
        if backgroundSize != 0:
            mean = (sum / N) - backgndMean
            sumcorrected = sum - (N * backgndMean)
            result = [sum, mean, N]
        if backgroundSize == 0:
            mean = (sum / N)
            result = [sum, mean, N]
        if output is True:
            print ('Sum:', result[0], 'Mean:', result[1],
                   'Number of Pixel:', result[2])
        if newAnnotation is False:
            fileout = open('apertures.ann', 'a')
        if newAnnotation is True:
            fileout = open('apertures.ann', 'w')
        if annotation is True:
            position = astFunc.equatorial_to_degrees(position)
            apertureString = ('CROSS ' +
                              str(position[0]) + ' ' +
                              str(position[1]) + ' 0.0008 0.001 \n'
                              'circle  ' +
                              str(position[0]) + ' ' +
                              str(position[1]) + ' ' +
                              str(apertureSize * (1. / 60 / 60) / 2) + ' \n'
                              'text ' +
                              str(position[0] + (1. / 60 / 60)) + ' ' +
                              str(position[1] + (1. / 60 / 60)) + '\n')
            fileout.write(apertureString)
            fileout.close()
        return result

    def sky2pix(self, coordinate, origin=0):
        r"""
        Calculates the Pixel corresponding to a given coordinate.

        Parameters
        -----------
        coordinate : list
            Either ['RA','DEC'], e.g. ['1:34:7.00', '+30:47:52.00'] in
            equatorial coordinates or [RA, DEC] in GRAD.
        origin : int
           ``0`` or ``1``; this steers how the first pixel is counted
           ``0`` is for usage with python as it starts to count from zero.
           ``1`` is the fits standart.

        Returns
        --------
        pixel : List
            [x, y]; the pixel coordinates of the map.
        """
        # Equatorial_to_degrees is only possible to execute if coordinate is in
        # the correct form for RA DEC.
        try:
            coordinate = astFunc.equatorial_to_degrees(coordinate)
            print coordinate
            print [[coordinate[0]],[coordinate[1]]]
        except:
            pass
        coordinate = np.array([[coordinate[0], coordinate[1]]], np.float_)
        pixel = self.wcs.wcs_sky2pix(coordinate, origin)[0]
        # The pixels are floats. To get integers we get the floor value
        # of the floats
        pixel = [int(math.floor(float(pixel[0]))),
                 int(math.floor(float(pixel[1])))]
        return pixel

    def pix2sky(self, pixel, degrees_or_equatorial='degrees'):
        r"""
        Calculates the Coordinates of a given Pixel.

        Parameters
        -----------
        pixel : list
            Pixel of the map; [x, y].
        degrees_or_equatorial : string
            Either ``"degrees"`` or ``"equatorial"``. Choosing the
            Format of the coordintates to be returnes.
            Defaults to ``"degrees"``.

        Returns
        --------
        coordinate : list
            The coordinates corresponding to pixel. Either in Degrees or in
            Equatorial coordinates, depending on the parameter
            *degrees_or_equatorial*.
        """
        # equatorial_to_degrees is only possible to execute if coordinate
        # is in RA DEC
        coordinate = self.wcs.wcs_pix2sky([[float(pixel[0]),
                                            float(pixel[1])]], 1)[0]
        if degrees_or_equatorial == 'equatorial':
            coordinate = astFunc.degrees_to_equatorial(coordinate)
        return coordinate

    def gauss_factor(self, beamConv, beamOrig=None, dx1=None, dy1=None):
        r"""
        Caluclates the scaling factor to be applied after convolving
        a map in Jy/beam with a gaussian to get fluxes in Jy/beam again.

        This function is a copy of the FORTRAN gaufac function from the Miriad
        package, which determine the Gaussian parameters resulting from
        convolving two gaussians. This function yields the same result as
        the MIRIAD gaufac function.

        Parameters
        ----------
        beamConv : list
            A list of the [major axis, minor axis, position_angle]
            of the gaussion used for convolution.
        beamOrig :
            Same format as beamConv but giving the parameters of the original
            beam of the map. As Default the self.resolution list is used.
        dx1, dy1 : floats
            Being the pixel size in both dimensions of the map.
            By default the ``CDELT1`` and ``CDELT2`` keywords from the
            fits header are used.

        Returns
        -------
        fac :
            Factor for the output Units.
        amp :
            Amplitude of resultant gaussian.
        bmaj, bmin :
            Major and minor axes of resultant gaussian.
        bpa :
            Position angle of the resulting gaussian.
        """
        #include 'mirconst.h'
        # Define cosine and Sinus of the position Angles of the
        # Gaussians
        bmaj2, bmin2, bpa2 = beamConv
        bmaj2, bmin2, bpa2 = (bmaj2 * const.arcsecInGrad, bmin2 *
                              const.arcsecInGrad, bpa2 * const.arcsecInGrad)
        if beamOrig is None:
            bmaj1, bmin1, bpa1 = self.resolution
            bmaj1, bmin1, bpa1 = (bmaj1 * const.arcsecInGrad,
                                  bmin1 * const.arcsecInGrad,
                                  bpa1 * const.arcsecInGrad)
        if dx1 is None:
            dx1 = self.header['CDELT1']
        if dy1 is None:
            dy1 = self.header['CDELT1']
        cospa1 = math.cos(bpa1)
        cospa2 = math.cos(bpa2)
        sinpa1 = math.sin(bpa1)
        sinpa2 = math.sin(bpa2)
        alpha = ((bmaj1 * cospa1) ** 2
                 + (bmin1 * sinpa1) ** 2
                 + (bmaj2 * cospa2) ** 2
                 + (bmin2 * sinpa2) ** 2)
        beta = ((bmaj1 * sinpa1) ** 2
                + (bmin1 * cospa1) ** 2
                + (bmaj2 * sinpa2) ** 2
                + (bmin2 * cospa2) ** 2)
        gamma = (2 * ((bmin1 ** 2 - bmaj1 ** 2)
                      * sinpa1 * cospa1
                      + (bmin2 ** 2 - bmaj2 ** 2)
                      * sinpa2 * cospa2))
        s = alpha + beta
        t = math.sqrt((alpha - beta) ** 2 + gamma ** 2)
        bmaj = math.sqrt(0.5 * (s + t))
        bmin = math.sqrt(0.5 * (s - t))
        if (abs(gamma) + abs(alpha - beta)) == 0:
            bpa = 0.0
        else:
            bpa = 0.5 * atan2(-1 * gamma, alpha - beta)
        #print math.pi/ (4.0*math.log(2.0))
        amp = (math.pi / (4.0 * math.log(2.0)) * bmaj1 * bmin1 * bmaj2 * bmin2
               / math.sqrt(alpha * beta - 0.25 * gamma * gamma))
        fac = ((math.sqrt(dx1 ** 2) * math.sqrt(dy1 ** 2))) / amp
        #print fac
        return fac, amp, bmaj * 60 * 60, bmin * 60 * 60, bpa

    def change_unit(self, final_unit, frequency=None, debug=False):
        r"""
        Changes the unit of a map in an automated way.

        Parameters
        ----------
        final_unit : string
            The unit to change the map to. Possible are:
                1. Jy/beam: ``"JyB"`` ``"JyBeam"``
                2. Jy/pixel: ``"JyP"``, ``"JyPix"``, ``"JyPixel"``
                3. MJy/sterad: ``"MJyPsr"``, ``"MJPSR"``, ``"MJy/sr"``
                4. Main Beam Temperature: ``"Tmb"``, ``"T"``, ``"Kkms"``
                5. erg/s/pixel: ``"ergs"`` ``"ERGPSECPPIX"``,
                                ``"ERGPSECPPIXEL"``, ``"ERG-S-1-Pixel"``,
                                ``"ERG-S-1"``
                6. erg/s/beam: ``"ERGPSECPBEAM"``
                7. erg/s/sterad ``"ERGPERSTER"``
        frequency : float
            Can be used if self.frequency is NaN. The frequency (in GHz) is
            needed for conversions between temperature and Jansky/Erg scale.
            Other conversions don't need it.

        Notes
        -----

        .. warning::

            This function is still in development and not all
            conversions may work properly.
        """
        # Definition of the unit nomenclatures.
        # TODO!! This should be Uset input.
        _jansky_beam_names = ['JYB', 'JyBeam']
        _jansky_pixel_names = ['JYPIX', 'JYP', 'JYPIXEL']
        _tmb_names = ['TMB', 'T', 'KKMS']
        _MJy_per_sterad_names = ['MJPSR', 'MJYPSR', 'MJY/SR']
        _erg_sec_pixel_names = ['ERGPSECPPIX', 'ERGPSECPPIXEL',
                                'ERG-S-1-Pixel', 'ERG-S-1']
        _erg_sec_beam_names = ['ERGPSECPBEAM']
        _erg_sec_sterad_names = ['ERGPERSTER']
        _known_units = (_jansky_beam_names + _jansky_pixel_names +
                        _tmb_names + _MJy_per_sterad_names)

        # If original and final unit are the same no conversion is needed.
        if final_unit.upper() == self.fluxUnit.upper():
                self.update_file()

        # Jy/Pixel to Jy/B.
        def JyPix_JyBeam(invert=False):
            '''
            This function converts between Jy/Pixel and Jy/Beam.

            Parameters
            ----------

            invert : logic
                If ``"True"`` the conversion is done from ``Jy/Beam``
                to `` Jy/Pixel``. If ``"False"`` vice-versa.
            '''
            if not invert:
                factor = self.beamSizeSterad / self.pixelSizeSterad
                self.data = self.data * factor
                self.fluxUnit = 'JyB'
                self.header.update('BUNIT', 'Jy/Beam')
            if invert:
                factor = self.pixelSize / self.beamSize
                self.data = self.data * factor
                self.fluxUnit = 'JyPix'
                self.header.update('BUNIT', 'Jy/Pixel')

        def MJyPSr_JyBeam(invert=False):
            if not invert:
                # Working
                # fwhm 24micron - 6", 70micron - 18", 160micron - 40"
                # MJy/sterad = 10^6 /(4.25*10^10/(1.133*fwhm[arcsec]^2)) J/beam
                # with 4.25 = 1/(4.848e-6)^2 -> 1'' = 4.848e-6 rad
                # - > Jy/beam = 2.6629e-5 * fwhm[arcsec]^2 MJy/sr
                factor = (2.6629e-5
                          * self.resolution[0]
                          * self.resolution[1])
                self.data = self.data * factor
                self.fluxUnit = 'JyB'
                self.header.update('BUNIT', 'Jy/Beam')
            if invert:
                # Working
                # fwhm 24micron - 6", 70micron - 18", 160micron - 40"
                # MJy/sterad = 10^6 /(4.25*10^10/(1.133*fwhm[arcsec]^2)) J/beam
                # with 4.25 = 1/(4.848e-6)^2 -> 1'' = 4.848e-6 rad
                # - > Jy/beam = 2.6629e-5 * fwhm[arcsec]^2 MJy/sr
                factor = (2.6629e-5
                          * self.resolution[0]
                          * self.resolution[1])
                factor = 1 / factor
                self.data = self.data * factor
                self.fluxUnit = 'MJpSr'
                self.header.update('BUNIT', 'MJy/sterad')

        def Tmb_JyBeam(invert=False):
            if not invert:
                if self.frequency is not np.nan:
                    self.data = self.data * self.flux_conversion()
                    self.fluxUnit = 'JyB'
                    self.header.update('BUNIT', 'Jy/Beam')
                elif frequency is not None:
                    factor = self.flux_conversion(frequency=frequency)
                    self.data = self.data * factor
                    self.fluxUnit = 'JyB'
                    self.header.update('BUNIT', 'Jy/Beam')
                else:
                    print ('Wavelength needed to convert between Jy and K.'
                           '\nCannot continue -> Exit!')
                    sys.exit()
            if invert:
                if self.frequency is not np.nan:
                    factor = self.flux_conversion()
                    self.data = self.data * factor
                    self.fluxUnit = 'Tmb'
                    self.header.update('BUNIT', 'K km/s')
                elif frequency is not None:
                    factor = self.flux_conversion(frequency=frequency)
                    self.data = self.data * factor
                    self.fluxUnit = 'Tmb'
                    self.header.update('BUNIT', 'K km/s')
                else:
                    print ('Wavelength needed to convert between Jy and K.'
                           '\nCannot continue -> Exit!')
                    sys.exit()

        def ErgSecBeam_JyBeam(invert=False):
            if not invert:
                factor = units.JyBToErgsB(None, self.distance, self.wavelength,
                                          invert=True, map_use=True)
                self.data = self.data * factor
                self.fluxUnit = 'JyB'
                self.header.update('BUNIT', 'Jy/beam')
            if invert:
                factor = units.JyBToErgsB(None, self.distance, self.wavelength,
                                          invert=False, map_use=True)
                self.data = self.data * factor
                self.fluxUnit = 'ErgPSecPBeam'
                self.header.update('BUNIT', 'ergs s^-1 beam^-1')

        def ErgSecBeam_ErgSecPixel(invert=False):
            if not invert:
                factor = self.beamSizeSterad / self.pixelSizeSterad
                self.data = self.data * factor
                self.fluxUnit = 'ErgPSecPPixel'
                self.header.update('BUNIT', 'ergs s^-1 pixel^-1')
            if invert:
                factor = self.pixelSizeSterad / self.beamSizeSterad
                self.data = self.data * factor
                self.fluxUnit = 'ErgPSecPBeam'
                self.header.update('BUNIT', 'ergs s^-1 beam^-1')

        # Conversions
        # Jy/pixel <-> Jy/beam
        if final_unit.upper() in _jansky_beam_names:
            if self.fluxUnit.upper() in _jansky_pixel_names:
                JyPix_JyBeam()
                self.update_file()

        if final_unit.upper() in _jansky_pixel_names:
            if self.fluxUnit.upper() in _jansky_beam_names:
                JyPix_JyBeam(invert=True)
                self.update_file()

        # Jy/beam <-> MJy/sterad.
        if final_unit.upper() in _jansky_beam_names:
            if self.fluxUnit.upper() in _MJy_per_sterad_names:
                MJyPSr_JyBeam()
                self.update_file()
        if final_unit.upper() in _MJy_per_sterad_names:
            if self.fluxUnit.upper() in _jansky_beam_names:
                MJyPSr_JyBeam(invert=True)
                self.update_file()

        # Jy/pixel <-> MJy/sterad.
        if final_unit.upper() in _jansky_pixel_names:
            if self.fluxUnit.upper() in _MJy_per_sterad_names:
                MJyPSr_JyBeam()
                JyPix_JyBeam(invert=True)
                self.update_file()
        if final_unit.upper() in _MJy_per_sterad_names:
            if self.fluxUnit.upper() in _jansky_pixel_names:
                JyPix_JyBeam()
                MJyPSr_JyBeam(invert=True)
                self.update_file()

        # MJy/sterad <-> Tmb
        if final_unit.upper() in _tmb_names:
            if self.fluxUnit.upper() in _MJy_per_sterad_names:
                MJyPSr_JyBeam()
                Tmb_JyBeam(invert=True)
                self.update_file()
        if final_unit.upper() in _MJy_per_sterad_names:
            if self.fluxUnit.upper() in _tmb_names:
                Tmb_JyBeam()
                MJyPSr_JyBeam(inver=True)
                self.update_file()

        # Jy/Beam <-> Tmb (Main beam temperature).
        if final_unit.upper() in _tmb_names:
            if self.fluxUnit.upper() in _jansky_beam_names:
                Tmb_JyBeam(invert=True)
                self.update_file()
        if final_unit.upper() in _jansky_beam_names:
            if self.fluxUnit.upper() in _tmb_names:
                Tmb_JyBeam()
                self.update_file()

        # Tmb <-> Jy/pixel
        if final_unit.upper() in _tmb_names:
            if self.fluxUnit.upper() in _jansky_pixel_names:
                JyPix_JyBeam()
                Tmb_JyBeam(invert=True)
                self.update_file()
        if final_unit.upper() in _jansky_pixel_names:
            if self.fluxUnit.upper() in _tmb_names:
                Tmb_JyBeam()
                JyPix_JyBeam(invert=True)
                self.update_file()

        # Jy/beam <-> Erg/s/beam
        if final_unit.upper() in _erg_sec_beam_names:
            if self.fluxUnit.upper() in _jansky_beam_names:
                ErgSecBeam_JyBeam(invert=True)
                self.update_file()
        if final_unit.upper() in _jansky_beam_names:
            if self.fluxUnit.upper() in _erg_sec_beam_names:
                ErgSecBeam_JyBeam()
                self.update_file()

       # Jy/beam <-> Erg/s/pixel
        if final_unit.upper() in _erg_sec_pixel_names:
            if self.fluxUnit.upper() in _jansky_beam_names:
                ErgSecBeam_JyBeam(invert=True)
                ErgSecBeam_ErgSecPixel()
                self.update_file()
        if final_unit.upper() in _jansky_beam_names:
            if self.fluxUnit.upper() in _erg_sec_pixel_names:
                ErgSecBeam_ErgSecPixel(invert=True)
                ErgSecBeam_JyBeam()
                self.update_file()

       # Jy/pixel <-> Erg/s/pixel
        if final_unit.upper() in _erg_sec_pixel_names:
            if self.fluxUnit.upper() in _jansky_pixel_names:
                JyPix_JyBeam()
                ErgSecBeam_JyBeam(invert=True)
                ErgSecBeam_ErgSecPixel()
                self.update_file()
        if final_unit.upper() in _jansky_pixel_names:
            if self.fluxUnit.upper() in _erg_sec_pixel_names:
                ErgSecBeam_ErgSecPixel(invert=True)
                ErgSecBeam_JyBeam()
                JyPix_JyBeam(invert=True)
                self.update_file()

       # Jy/pixel <-> Erg/s/beam
        if final_unit.upper() in _erg_sec_beam_names:
            if self.fluxUnit.upper() in _jansky_pixel_names:
                JyPix_JyBeam()
                ErgSecBeam_JyBeam(invert=True)
                self.update_file()
        if final_unit.upper() in _jansky_pixel_names:
            if self.fluxUnit.upper() in _erg_sec_beam_names:
                ErgSecBeam_JyBeam()
                JyPix_JyBeam(invert=True)
                self.update_file()

        # MJy/sterad <-> Erg/s/pixel
        if final_unit.upper() in _erg_sec_pixel_names:
            if self.fluxUnit.upper() in _MJy_per_sterad_names:
                MJyPSr_JyBeam()
                ErgSecBeam_JyBeam(invert=True)
                ErgSecBeam_ErgSecPixel()
                self.update_file()
        if final_unit.upper() in _MJy_per_sterad_names:
            if self.fluxUnit.upper() in _erg_sec_pixel_names:
                ErgSecBeam_ErgSecPixel(invert=True)
                ErgSecBeam_JyBeam()
                MJyPSr_JyBeam(invert=True)
                self.update_file()

        # MJy/sterad <-> Erg/s/beam
        if final_unit.upper() in _erg_sec_beam_names:
            if self.fluxUnit.upper() in _MJy_per_sterad_names:
                MJyPSr_JyBeam()
                ErgSecBeam_JyBeam(invert=True)
                self.update_file()
        if final_unit.upper() in _MJy_per_sterad_names:
            if self.fluxUnit.upper() in _erg_sec_beam_names:
                ErgSecBeam_JyBeam()
                MJyPSr_JyBeam(invert=True)
                self.update_file()

        # MJy/sterad <-> Erg/s/beam
        if final_unit.upper() in _erg_sec_pixel_names:
            if self.fluxUnit.upper() in _erg_sec_beam_names:
                ErgSecBeam_ErgSecPixel()
                self.update_file()
        if final_unit.upper() in _erg_sec_beam_names:
            if self.fluxUnit.upper() in _erg_sec_pixel_names:
                ErgSecBeam_ErgSecPixel(invert=True)
                self.update_file()

       # Tmb <-> Erg/s/pixel
        if final_unit.upper() in _erg_sec_pixel_names:
            if self.fluxUnit.upper() in _tmb_names:
                Tmb_JyBeam()
                ErgSecBeam_JyBeam(invert=True)
                ErgSecBeam_ErgSecPixel()
                self.update_file()
        if final_unit.upper() in _tmb_names:
            if self.fluxUnit.upper() in _erg_sec_pixel_names:
                ErgSecBeam_ErgSecPixel(invert=True)
                ErgSecBeam_JyBeam()
                Tmb_JyBeam(invert=True)
                self.update_file()

       # Tmb <-> Erg/s/beam
        if final_unit.upper() in _erg_sec_beam_names:
            if self.fluxUnit.upper() in _tmb_names:
                Tmb_JyBeam()
                ErgSecBeam_JyBeam(invert=True)
                self.update_file()
        if final_unit.upper() in _tmb_names:
            if self.fluxUnit.upper() in _erg_sec_beam_names:
                ErgSecBeam_JyBeam()
                Tmb_JyBeam(invert=True)
                self.update_file()

        # If the original flux unit is not known.
        elif self.fluxUnit.upper() not in _known_units:
            print ('Flux unit not known to astrolyze - '
                   'Can not continue')
            if debug is True:
                pass
            elif debug is False:
                sys.exit()
        return FitsMap(self.returnName())

    def toGildas(self, prefix=None):
        r"""
        Changes the current map to the Gildas Format.

        The function takes changes to the map_name variables
        made outside of functions into account via
        :py:func:`astrolyze.maps.main.Map.returnName` into account.

        Parameters
        ----------
        prefix: string or None
            Path to location where the new gildas file will be stored.
            The default is None which defaults to the current self.prefix.

        Examples
        --------

        To continue working with the gildas map use:

        >>> map = map.toGildas()

        To only store the current map in the gildas format and go on
        working with the fits file use:

        >>> map.toGildas()

        Here map is an Instance of the FitsMap class.
        """
        prefix = prefix or self.prefix
        gildas_name = self.returnName(prefix=prefix, dataFormat='gdf')
        _conv_file = open('temp.greg', 'w')
        _conv_file.write('fits ' + self.map_name + ' to ' + gildas_name + '\n'
                         'exit\n')
        _conv_file.close()
        os.system('greg -nw @temp.greg')
        os.system('rm temp.greg')
        return gildas.GildasMap(gildas_name)

    def toMiriad(self, prefix=None):
        r"""
        Changes the current map to the Miriad Format.

        The function takes changes to the map_name variables made outside of
        functions into account via :py:func:`maps.main.Map.returnName` into
        account.

        Parameters
        ----------
        prefix : string or None
            Path to location where the new gildas file will be stored.
            The default is None which defaults to the current self.prefix.

        Examples
        --------

        This function works like
        :py:func:`maps.mapClassFits.FitsMap.toGildas` and the same
        Examples apply.

        """
        _prefix = prefix or self.prefix
        _fileout = open('miriad.out', 'a')
        _miriad_name = self.returnName(prefix=_prefix, dataFormat=None)
        os.system('rm -rf ' + _miriad_name)
        _fileout.write('fits in=' + str(self.map_name) + ' out=' +
                       _miriad_name + ' op=xyin\n')
        os.system('fits in=' + str(self.map_name) + ' out=' +
                  _miriad_name + ' op=xyin')
        os.system('rm miriad.out')
        return miriad.MiriadMap(_miriad_name)

    def toFits(self):
        if self.map_name != self.returnName():
            os.system('cp -rf ' + self.map_name + ' ' + self.returnName())
            return fits.FitsMap(self.returnName())
        else:
            return self
