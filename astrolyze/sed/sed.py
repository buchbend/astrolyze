# Copyright (C) 2012, Christof Buchbender
# BSD Licencse
import sys
import os
import numpy as np
import pyfits
import matplotlib.pyplot as plt

import astrolyze.maps.main as main
import astrolyze.maps.fits as fits
import astrolyze.maps.gildas as gildas
import astrolyze.maps.miriad as miriad
import astrolyze.maps.stack as stack

import astrolyze.maps.tools as mtools
import astrolyze.functions.astro_functions as astro_functions


class SedStack(stack.Stack):
    r"""
    Reads in the SEDs from the maps stored under the input folder at given
    coordinates and creates a stack of Sed objects.
    """
    def __init__(self, folder, data_format='.fits', coordinates=False,
                 flux_acquisition='aperture', aperture=120, annotation=False,
                 full_map=False, output_folder=None, number_components=2):
        # Read in the stack from the input folder.
        stack.Stack.__init__(self, folder=folder, data_format=data_format)
        # Pass the input parameters to the variables
        self.flux_acquisition = flux_acquisition
        self.aperture = aperture
        self.number_components = number_components
        if self.flux_acquisition == 'pixel':
            if max(self.resolutions) != min(self.resolutions):
                print 'Maps do not have the same resolutions'
                sys.exit()
        self.annotation = annotation
        # check if the coordinates are given as
        if not coordinates and full_map and output_folder:
            self.get_map_seds(output_folder)
        if coordinates and not full_map:
            if type(coordinates) is str:
                self.load_coordinates(coordinates)
            if type(coordinates) is list:
                self.coordinates = coordinates
            self.get_seds()

    def info(self):
        for sed in self.sed_stack:
            sed.info()

    def load_coordinates(self, input_file):
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
        self.source_names: list
            Format::
            [[source_name_1, RA_1, DEC_1] , ... , [source_name_N, RA_N, DEC_N]]

        Examples
        --------

        The format of the coordinates given in the file must be in
        Equatorial:
        >>> equatorial_coordinates = ['02:23:21', '-02:23:21']
        """
        lines = open(input_file).readlines()
        self.coordinates = []
        self.source_names = []
        for line in lines:
            line_items = line.split()
            self.source_names += [line_items[0]]
            self.coordinates += [[line_items[1], line_items[2]]]
        if len(self.source_names) != len(self.coordinates):
            # source_names and coordinates have to be the same
            # length.
            # ADD RAISE Exception correctly
            raise ValueError("Something went wrong reding the cooridnates.")
        self.number_of_seds = len(self.coordinates)

    def get_seds(self):
        r""" Creates a stack of SEDs from the stack of maps loaded from the
        input folder if particular coordinates are given.
        """
        self.sed_stack = []
        for x, coordinate in enumerate(self.coordinates):
            flux_array = [[], [], []]
            names = []
            for maps in self.stack:
                if self.flux_acquisition == 'aperture':
                    flux = maps.read_aperture(coordinate,
                                         apertureSize=self.aperture,
                                         annotation=self.annotation)
                    flux_array[0] += [maps.frequency/1e9]
                    flux_array[1] += [flux[0]]
                    flux_array[2] += [flux[0] * maps.calibrationError]
                if self.flux_acquisition == 'pixel':
                    flux = maps.read_flux(coordinate)
                    flux_array[0] += [maps.frequency/1e9]
                    flux_array[1] += [flux]
                    flux_array[2] += [flux * maps.calibrationError]
                names += self.source_names[x]
            flux_array = np.asarray(flux_array)
            self.sed_stack += [Sed(self.source_names[x], coordinate,
                                   flux_array, self.number_components)]

    def get_map_seds(self, folder):
        r""" This functions fits the SED at every pixel of the input maps.
        Note that the maps have to be exactly the same size for this function to
        work. This can be achieved with e.g.::

          from astrolyze import *
          stack = Stack('some_input_folder')
          folder = 'some_output_folder'
          template = 'Path_to_a _template_map'
          stack = stack.unify_dimension(template, folder)

        Parameters
        ----------

        folder: string
            The path to the folder where the temperature, mass, beta and chisq
            maps are created. If the folder does not exist is will be created.

        Notes
        -----

        Depending on the number of pixels in the input image this function may
        take a good while to finish.
        """
        def progress(x):
            out = '%1.1f %% Done.' % x  # The output
            bs = '\b' * 1000            # The backspace
            print bs,
            print out,
        source_name = self.stack[0].source
        # Assure that the data array has only 2 dimensions. Sometimes there are
        # higher unused dimensions stores int he fits files.
        while len(self.stack[0].data) == 1:
            self.stack[0].data = self.stack[0].data[0]
        # Initialization of the arrays that store the fitted values at every
        # points of the maps. We make a copy of an original data array.
        self.temperature_maps = []
        self.mass_maps = []
        for i in range(self.number_components):
            self.temperature_maps += [self.stack[0].data.copy()]
            self.mass_maps += [self.stack[0].data.copy()]
        self.beta_map = self.stack[0].data.copy()
        self.chisq_map = self.stack[0].data.copy()
        # Now we walk over all pixels of the maps create the seds at every
        # point, fit it and write the results to the corresponding data_arrays.
        print 'Fitting the SEDs:'
        for i in range(len(self.stack[0].data)):
            percentage = 100./len(self.stack[0].data) * i
            progress(percentage)
            for j in range(len(self.stack[0].data[0])):
                flux_array = [[], [], []]
                for k, map_ in enumerate(self.stack):
                    while len(map_.data) == 1:
                        map_.data = map_.data[0]
                    flux_array[0] += [map_.frequency/1e9]
                    flux_array[1] += [map_.data[i][j]]
                    flux_array[2] += [map_.data[i][j] *
                                      map_.calibrationError]
                flux_array = np.asarray(flux_array)
                coordinate = [i, j]
                sed = Sed(source_name, coordinate, flux_array,
                    self.number_components, init_fit=True)
                for x, temp in enumerate(sed.fit_temperatures):
                    self.temperature_maps[x][i][j] = temp
                for x, mass in enumerate(sed.fit_masses):
                    self.mass_maps[x][i][j] = mass
                self.beta_map[i][j] = sed.fit_beta[0]
                self.chisq_map[i][j] = sed.fit_chisq
        # Finally we create fits maps from the stored results and save them in
        # the given output folder. It the latter does not exist we create it.
        print 'Creating Maps ...'
        if folder is not None and '/' not in folder[-1]:
            folder = folder + '/'
        if not os.path.isdir(folder):
            os.system('mkdir ' + folder)
        header = self.stack[0].header
        for x, item in enumerate(self.temperature_maps):
            header['BUNIT'] = 'Kelvin'
            header['DATAMIN'] = item.min()
            header['DATAMAX'] = item.max()
            pyfits.writeto('{0}{1}_SED_Temp{2}_Kelvin_{3}'
                           '.fits'.format(folder, str(self.stack[0].source),
                                       str(x+1),
                                       str(self.stack[0].resolutionToString())),
                           item, header)
        for x, item in enumerate(self.mass_maps):
            header['BUNIT'] = 'Msun'
            header['DATAMIN'] = item.min()
            header['DATAMAX'] = item.max()
            pyfits.writeto('{0}{1}_SED_Mass{2}_Msun_{3}'
                           '.fits'.format(folder, str(self.stack[0].source),
                                       str(x+1),
                                       str(self.stack[0].resolutionToString())),
                           item, self.stack[0].header)
        header['BUNIT'] = ''
        header['DATAMIN'] = self.beta_map.min()
        header['DATAMAX'] = self.beta_map.max()
        pyfits.writeto('{0}{1}_SED_Beta_None_{2}'
                       '.fits'.format(folder, str(self.stack[0].source),
                                      str(self.stack[0].resolutionToString())),
                       self.beta_map, self.stack[0].header)
        header['DATAMIN'] = self.chisq_map.min()
        header['DATAMAX'] = self.chisq_map.max()
        pyfits.writeto('{0}{1}_SED_Chisq_None_{2}'
                       '.fits'.format(folder, str(self.stack[0].source),
                                      str(self.stack[0].resolutionToString())),
                       self.chisq_map, self.stack[0].header)
        print 'Finished!'


class Sed:
    r""" This class handles a single SED. Basically it is able to fit,
    and plot the sed.

    Parameters
    ----------

    source_name: string
        The name of the source to which the SED corresponds to.
    coordinate: list
        The coordinate of the source. [RA, DEC]
    flux_array: list
        The array that is created by SedStack with the entries of wavelength,
        flux, and error.
    init_fit: logic
        Steers whether the SED is fitted already during creation.
    """
    def __init__(self, source_name, coordinate, flux_array,
                 number_components=2, init_fit=True):
        self.source_name = source_name
        self.coordinate = coordinate
        self.number_components = number_components
        self.flux_array = flux_array
        self.fit_temperatures = None
        self.fit_masses = None
        self.fit_beta = None
        self.fit_done = False
        self.set_defaults()
        if init_fit:
            self.grey_body_fit()

    def info(self):
        print 'Source Name: ' + str(self.source_name)
        print 'Coordinate: ' + str(self.coordinate)
        print 'flux_array: ' + str(self.flux_array)
        if self.fit_done:
            print 'Temperature Fit: ' + str(self.fit_temperatures)
            print 'Mass Fit: ' + str(self.fit_masses)
            print 'Beta Fit: ' + str(self.fit_beta)
            print 'Chisq: ' + str(self.fit_chisq)
        if not self.fit_done:
            print 'SED not fitted yet'

    def set_defaults(self):
        # Set a default guess for the input parameter
        # needed to run a greybody Fit.
        init_guess_temperature = [20, 50, 100, 150, 200]
        init_guess_masses = [1e5, 1e2, 1e1 ,1e-1, 1e-2]
        init_beta_guess = [2.]
        self.temperature_guess = []
        self.mass_guess = []
        for x in range(self.number_components):
            self.temperature_guess += [init_guess_temperature[x]]
            self.mass_guess += [init_guess_masses[x]]
        self.beta_guess = [2.]
        self.p1 = [self.temperature_guess, self.mass_guess, self.beta_guess]
        self.p2 = None
        # Setting up default input choices for the
        # grey_body_fit function from astroFunctions.
        # please check this function for more details.
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

        ..  :py:func:`astrolyze.functions.astro_functions.grey_body_fit`
        """
        try:
            (self.p2, 
             self.fit_chisq) = astro_functions.grey_body_fit(
                                           data=self.flux_array,
                                           start_parameter=self.p1,
                                           nu_or_lambda=self.nu_or_lambda,
                                           fit_beta=self.fit_beta,
                                           fix_temperature=self.fix_temperature,
                                           kappa=self.kappa,
                                           residuals=self.residuals)
        except:
            raise ValueError('Data could not be fitted!')
        self.fit_temperatures = self.p2[0]
        self.fit_masses = self.p2[1]
        self.fit_beta = self.p2[2]
        self.fit_done = True
        # print 'Fit-Done'

    def plot_sed(self, axes=plt.gca(), nu_or_lambda='nu', color='black',
                 linewidth=0.5, xRange='normal'):
        '''Plot a multi component greybody model.

        Parameters
        ----------

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
            model,grey = astro_functions.multi_component_grey_body(self.p2, x,
                                                                   'nu',
                                                                   self.kappa)
        if nu_or_lambda=='lambda':
            model,grey = astro_functions.multi_component_grey_body(self.p2, x,
                                                                   'nu',
                                                                   self.kappa)
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
        plt.loglog(x,model,ls='-', color=color, label='_nolegend_', lw=0.5,
                   marker='')
        linestyles = [':','-.','-']
        j=0
        for i in grey:
            plt.loglog(x, i, color=color, ls=linestyles[j], lw=0.5, marker='')
            j+=1
