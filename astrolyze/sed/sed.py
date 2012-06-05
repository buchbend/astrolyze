# Copyright (C) 2012, Christof Buchbender
# BSD Licencse
import sys
import numpy as np
import matplotlib.pyplot as plt

import astrolyze.maptools as mt
import astrolyze.functions.astro_functions as af

from astrolyze.maps import *

class SedStack:
    r"""
    Reads in the SEDs from the maps stored under the input folder at given
    coordinates and creates a stack of Sed objects.
    """
    def __init__(self, folder, data_format='.fits', filein=None,
                 coordinates=None, flux_aquisition='aperture', aperture=120,
                 annotation=False):
        #Some Defaults
        self.flux_aquisition = flux_aquisition
        self.aperture = aperture
        self.annotation = annotation
        self.load_file_list(folder, data_format)
        if filein is None and coordinates is None: 
            pass
        if filein is None:
            self.load_coordinates(filein)
        if coordinates is None:
            self.coordinates = cooridnates
        self.get_seds()

    def load_file_list(self, folder, data_format):
        r"""
        Intitializes the variable self.file_list containing a list of the maps
        that are to be used to create the SED.

        Parameters
        ----------

        folder: string
            May contain sub-folders that are all taken into account.
        data_format: string
            Supported ``'.fits'`` (Default) and ``'.gdf'``. Maybe 
            extended later.

        Returns
        -------

        self.file_list: list
            Paths to the file in potential subfolders in under the folder.
        self.maps: list
            A list with FitsMap objects. If Gildas maps are loaded they are
            converted into fits maps.

        See Also
        --------

        maps.main, maps.fits, maps.gildas
        """
        self.file_list = mt.get_list(folder, data_format)
        self.maps = []
        # loading the maps with the mapclass
        for item in self.file_list:
            if '.fits' in item:
                map = FitsMap(item)
                self.maps += [map]
            if '.gdf' in item:
                map = GildasMap(item)
                map = map.toFits()
                self.maps += [map]

    def load_coordinates(self, filein):
        r"""
        Loads the coordinated where the SEDs should be evaluated
        from either a file or a list. Both are not possible.

        Parameters
        ----------

        filein: string
            Path to file that cotains the coordinates format has to be:

            source_name RA DEC

            RA, DEC has to be for epoch J2000 in Equatorial coordinates, 
            see below for examples of the syntax.

        Returns
        -------

        self.coordinates: list
            Format::
            [[source_name_1, RA_1, DEC_1] , ... , [source_name_N, RA_N, DEC_N]]

        Examples
        --------

        The format of the coordinates given in the file must be in
        Equatorial:
        >>> equatorial_coordinates = ['02:23:21', '-02:23:21']
        """
        lines = open(filein).readlines()
        self.coordinates = []
        self.source_names = []
        for line in lines:
            line_items = line.split()
            self.source_names += [line_items[0]]
            self.coordinates += [[line_items[1], line_items[2]]]
        if len(self.source_names) != len(self.coordinates):
            # source_names and coordinates have to be the same
            # lenght.
            # ADD RAISE Exception correctly
            self.number_of_seds = len(coordinates)

    def get_seds(self):
        self.sed_stack = []
        for x, coordinate in enumerate(self.coordinates):
            flux_array = [[], [], []]
            for map in self.maps:
                try:
                    if self.flux_aquisition_type == 'aperture':
                        flux = map.read_aperture(coordinate, 
                                         apertureSize=self.aperture, 
                                         annotation=self.annotation)
                    if self.flux_aquisition_type == 'pixel':
                        flux = map.read_flux(coordinate)
                    print  map.calibrationError
                    print map.frequency
                    flux_array[0] += [map.frequency/1e9]
                    flux_array[1] += [flux[0]]
                    flux_array[2] += [flux[0] * map.calibrationError]
                    names += self.source_names[x]
                except:
                    continue
            flux_array = np.asarray(flux_array)
            self.sed_stack += [Sed(coordinate, self.source_names[x], 
                                   flux_array)]


class Sed:
    r"""
    This class handles a single sed. Fitting plotting and so on. 
    It should contain the data, fitting, 
    """
    def __init__(self, source_name, coordinate, flux_array):
        self.source_name = source_name
        self.coordinate = coordinate
        self.flux_array = flux_array
        self.fit_done = False
        self.set_defaults()

    def set_defaults(self):
        # Set a default guess for the input parameter
        # needed to run a greybody Fit.
        self.temperature_guess = [20., 50.]
        self.mass_guess = [1e5, 1e2]
        self.beta_guess = [2.]
        self.p1 = [self.temperature_guess, self.mass_guess, self.beta_guess]
        self.p2 = None 
        # Setting up default input choices for the 
        # grey_body_fit function from astroFunctions.
        # please check there for their meaning.
        self.fit_beta = False
        self.kappa = 'Kruegel'
        self.fix_temperature = False
        self.rawChiSq = None
        self.residuals = False
        self.nu_or_lambda = 'nu'

    def grey_body_fit(self):
        r""""
        Fitting a multi componenet grey body to the input data in flux_array.

        See Also
        --------

..        :py:func:`astrolyze.functions.astro_functions.grey_body_fit`
        """
        try:
            self.p2, self.chisq = af.grey_body_fit(data=self.flux_array,
                                          start_parameter = self.p1, 
                                          nu_or_lambda=self.nu_or_lambda, 
                                          fit_beta=self.fit_beta, 
                                          fix_temperature=self.fix_temperature,
                                          kappa=self.kappa, 
                                          residuals=self.residuals)
            self.fit_temperatures = self.p2[0]
            self.fit_masses = self.p2[1]
            self.fit_beta = self.p2[2]
            self.fit_done = True
        except:
            pass
            #raise ValueError('Data could not be fitted!')
            #sys.exit()

    def plot_sed(self, axes=plt.gca(), nu_or_lambda='nu', color='black', linewidth=0.5, 
                 xRange='normal'):
        '''
        Plot a multi component greybody model.
 
        nu_or_lambda:
           plot against frequency ``'nu'`` or wavelenght ``'lambda'``
        kappa: 
            The kappa to use. ``'easy'`` or ``'Kruegel'``. Please refer
            to :py:func:`functions.astroFunctions.greyBody` for more
            information.
        xRange: PLEASE ADD DESCRIPTION
        linewidth: float
            The linewidth of the plotted lines. Default to 0.5.
        color: matplotlib conform color
            the color of the plotted lines. Default to ``'black'``.
        '''
        if not self.fit_done:
            raise ValueError('Could not plot the data. The data has not been '
                             'fitted yet.')
            sys.exit()
        if self.p2 ==  None:
            pass
        if xRange == 'LTIR':
        # Plot the SED in the range of the determination
            # of the L_TIR: 3-1100 micron 
            xmin =  3e-6# micron
            xmax =  4000e-6 # micron
            # conversion to frequency in GHz
            xmin = const.c/xmax/1e9
            xmax = const.c/xmin/1e9
            step = 0.1

        if xRange == 'normal':
            # arbitrary range definition 
            xmin = 1e-2
            xmax = 3e5	
            step = 0.5
        if type(xRange) == list:
            xmin = xRange[0]
            xmax = xRange[1]
            if len(xRange) < 3:
                step = 0.1
            else:
                step = xRange[2]
        x = np.arange(xmin,xmax,step)
        # multi_component_grey_body gives the summed 'model' and the components
        # grey'. 'grey' is a List 
        if nu_or_lambda == 'nu':
            model,grey = af.multi_component_grey_body(self.p2, x,'nu', self.kappa)
        if nu_or_lambda=='lambda':
            model,grey = af.multi_component_grey_body(self.p2, x,'nu', self.kappa)
            y=copy(x)
            modelLambda =copy(model)
            greyLambda = [] 
            for i in range(len(grey)):
                greyLambda += [copy(grey[i])]
            for i in range(len(x)):
                y[i]=(const.c/(x[len(x)-i-1]*1e9))/1e-6
            for i in range(len(model)):
                modelLambda[i]=model[len(model)-i-1]
                for j in range(len(greyLambda)):
                    greyLambda[j][i] = grey[j][len(grey[j])-i-1]
            x=y
            model =modelLambda
            grey = greyLambda
        plt.loglog(x,model,ls='-', color=color, label='_nolegend_', lw=0.5, marker='')
        linestyles = [':','-.','-']
        j=0
        for i in grey:
            plt.loglog(x, i, color=color, ls=linestyles[j], lw=0.5, marker='')
            j+=1

# Old Code


def getCoordinates(file):
    '''
    Returns a list of coordinates read from a file
    formated as follows:

        SourceName 1 RA1 DECN
        SourceName 2 RA1 DECN
        .
        .
        .
        SourceNameN RA DEC

    Input:  file -> i.e. File name
    Output: list with [[SourceName 1,...,SourceName N], [[RA 1,DEC 1],...,[RA
        N,DEC N]] ]
    '''
    filein = open(file).readlines()
    coordinates = [[], []]
    for i in filein:
        i = i.split()
        coordinates[0] += [i[0]]
        coordinates[1] += [[i[1], i[2]]]
    return coordinates

