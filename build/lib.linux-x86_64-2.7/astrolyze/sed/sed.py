# Copyright (C) 2012, Christof Buchbender
# BSD Licencse
import sys
import os
import numpy as np
import astropy.io.fits
import matplotlib.pyplot as plt

from generaltools import log_tools

import astrolyze.maps.main as main
import astrolyze.maps.fits as fits
import astrolyze.maps.gildas as gildas
import astrolyze.maps.miriad as miriad
import astrolyze.maps.stack as stack

import astrolyze.maps.tools as mtools
import astrolyze.functions.astro_functions as astro_functions
import astrolyze.functions.constants as const

USER = os.getenv("USER")

class SedStack(stack.Stack):
    r"""
    Reads in the SEDs from the maps stored under the input folder at given
    coordinates and creates a stack of Sed objects.

    Parameters
    ----------

    folder : str
        Input Folder
    data_format : str
        Filter on data format. Default `.fits`
    coordinates : list
        Coordinates where the SEDs should be extracted.
    full_map : bool
        If True temperature, mass and beta maps will be created by fitting the
        SEDs at every pixel.
    flux_acquisition : str
        `aperture` or `pixel` see explanation in
        :py:mod:`astrolyze.fits.FitsMap.flux_read`
    aperture : int
        Aperture to use for flux extraction if flux_acquisition == `aperture`
    annotation : bool
        If true a kvis annotation file is created with the positions where the
        SEDs, have been extracted.
    output_folder : str
        Where the data is stored.
    number_components : int
        Number of SED components to fit.
    mass_guess : list
        List with guesses for the initial dust mass. Has to have at
        least the same lenght as there are `number_components`
    temperature_guess : list
        List with guesses for the initial dust temperature. Has to have at
        least the same lenght as there are `number_components`
    beta_guess : list
        List with a guess for the beta value to be used has to be list with
        singel entry beta_guess=[2.].

    """
    def __init__(self, folder, data_format='.fits', coordinates=False,
                 flux_acquisition='pixel', aperture=None, annotation=False,
                 full_map=False, output_folder=None, number_components=2,
                 mass_guess=None, temperature_guess=None, beta_guess=None):
        # Read in the stack from the input folder.
        stack.Stack.__init__(self, folder=folder, data_format=data_format)
        # Pass the input parameters to the variables
        self.flux_acquisition = flux_acquisition
        self.aperture = aperture
        self.number_components = number_components
        self.mass_guess = mass_guess
        self.temperature_guess = temperature_guess
        self.beta_guess = beta_guess
        if self.flux_acquisition == 'pixel':
            if max(self.resolutions) != min(self.resolutions):
                self.log.error('Maps do not have the same resolutions')
                sys.exit()
        self.annotation = annotation
        # check if the coordinates are given as
        if coordinates and full_map:
            self.log.error("coordinates and full_map are exclusive")
            raise AttributeError
        if not coordinates and full_map and output_folder:
            self.get_map_seds(output_folder)
        if coordinates and not full_map:
            if type(coordinates) is str:
                self.load_coordinates(coordinates)
            if type(coordinates) is list:
                self.coordinates = coordinates
            self.get_seds()

    def __str__(self):
        string = "SedStack({})\nContains:\n".format(self.folder)
        for map_ in self.stack:
            string += "{}\n".format(map_)
        return string

    def __repr__(self):
        string = "SedStack({})".format(self.folder)
        return string

    def __len__(self):
        return len(self.stack)

    def __getitem__(self, position):
        return self.stack[position]

    def __iter__(self):
        return (i for i in self.stack)

    def info(self):
        for sed in self.sed_stack:
            sed.info()

    def load_coordinates(self, input_file):
        r"""
        Loads the coordinated where the SEDs should be evaluated
        from either a file or a list. Both are not possible.

        Parameters
        ----------

        filein : string
            Path to file that cotains the coordinates format has to be:

            source_name RA DEC

            RA, DEC has to be for epoch J2000 in Equatorial coordinates,
            see below for examples of the syntax.

        Returns
        -------

        self.coordinates : list
        self.source_names : list
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
                                   flux_array, self.number_components,
                                   mass_guess=self.mass_guess,
                                   temperature_guess=self.temperature_guess,
                                   beta_guess=self.beta_guess)]

    def get_map_seds(self, output_folder):
        r""" This functions fits the SED at every pixel of the input maps.

        Parameters
        ----------

        output_folder : string
            The path to the folder where the temperature, mass, beta and chisq
            maps are created. If the folder does not exist is will be created.

        Notes
        -----

        Depending on the number of pixels in the input image this function may
        take a good while to finish.

        Examples
        --------
        Note that the maps have to be exactly the same size for this function
        to work. This can be achieved with e.g.::

          from astrolyze import *
          stack = Stack('some_input_folder')
          output_folder = 'some_output_folder'
          template = 'Path_to_a _template_map'
          stack = stack.unify_dimension(template, folder)

        """
        def progress(x):
            out = '%1.1f %% Done.' % x  # The output
            bs = '\b' * 1000            # The backspace
            print bs,
            print out,
        source_name = self.stack[0].source
        # Assure that the data array has only 2 dimensions. Sometimes there are
        # higher unused dimensions stores in the fits files.
        while len(self.stack[0].data) == 1:
            self.stack[0].data = self.stack[0].data[0]
        # Initialization of the arrays that store the fitted values at every
        # point of the maps. We make a copy of an original data array.
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
                          self.number_components, init_fit=True,
                          temperature_guess=self.temperature_guess,
                          mass_guess=self.mass_guess,
                          beta_guess=self.beta_guess)
                for x, temp in enumerate(sed.fit_temperatures):
                    self.temperature_maps[x][i][j] = temp
                for x, mass in enumerate(sed.fit_masses):
                    self.mass_maps[x][i][j] = mass
                self.beta_map[i][j] = sed.fit_beta[0]
                self.chisq_map[i][j] = sed.fit_chisq
        # Finally we create fits maps from the stored results and save them in
        # the given output folder. It the latter does not exist we create it.
        print 'Creating Maps ...'
        if output_folder is not None and '/' not in output_folder[-1]:
            output_folder = output_folder + '/'
        if not os.path.isdir(output_folder):
            os.system('mkdir ' + output_folder)
        header = self.stack[0].header
        for x, item in enumerate(self.temperature_maps):
            header['BUNIT'] = 'Kelvin'
            header['DATAMIN'] = item.min()
            header['DATAMAX'] = item.max()
            astropy.io.fits.writeto('{0}{1}_SED_Temp{2}_Kelvin_{3}'
                           '.fits'.format(
                               output_folder, str(self.stack[0].source),
                               str(x+1),
                               str(self.stack[0].resolutionToString())),
                         item, header)
        for x, item in enumerate(self.mass_maps):
            header['BUNIT'] = 'Msun'
            header['DATAMIN'] = item.min()
            header['DATAMAX'] = item.max()
            header=self.stack[0].header
            astropy.io.fits.writeto('{0}{1}_SED_Mass{2}_Msun_{3}'
                           '.fits'.format(
                               output_folder, str(self.stack[0].source),
                               str(x+1),
                               str(self.stack[0].resolutionToString())),
                            item, header)
        header['BUNIT'] = ''
        header['DATAMIN'] = self.beta_map.min()
        header['DATAMAX'] = self.beta_map.max()
        data = self.beta_map
        header = header=self.stack[0].header
        astropy.io.fits.writeto(
            '{0}{1}_SED_Beta_None_{2}' \
            '.fits'.format(output_folder, str(self.stack[0].source),
                           str(self.stack[0].resolutionToString())),
            data, header
        )
        header['DATAMIN'] = self.chisq_map.min()
        header['DATAMAX'] = self.chisq_map.max()
        data = self.chisq_map
        header = header=self.stack[0].header
        astropy.io.fits.writeto(
            '{0}{1}_SED_Chisq_None_{2}' \
            '.fits'.format(output_folder, str(self.stack[0].source),
                           str(self.stack[0].resolutionToString())),
            data, header
        )
        print 'Finished!'


class Sed(object):
    r""" This class handles a single SED. Basically it is able to fit,
    and plot the sed.

    Parameters
    ----------

    source_name : string
        The name of the source to which the SED corresponds to.
    coordinate : list
        The coordinate of the source. [RA, DEC]
    flux_array : list
        The array that is created by SedStack with the entries of wavelength,
        flux, and error.
    init_fit : logic
        Controls whether the SED is fitted already during creation.
    number_components : int
        The number of greybody components to be fitted. Default: 2.
    temperature_guess : list
        List with the intial guesses for the temperature
    mass_guess : list
        List with the intial guesses for the masses
    beta_guess : list
        List with the intial guesses of the beta value
    """
    def __init__(self, source_name, coordinate, flux_array,
                 number_components=2, init_fit=True, temperature_guess=None,
                 mass_guess=None, beta_guess=None):
        self.log = log_tools.init_logger(
            directory="/home/{}/.astrolyze/".format(USER),
            name="astrolyze"
        )
        self.source_name = source_name
        self.coordinate = coordinate
        self.number_components = number_components
        self.temperature_guess = temperature_guess
        self.mass_guess = mass_guess
        self.beta_guess = beta_guess
        self.flux_array = flux_array
        self.fit_temperatures = None
        self.fit_masses = None
        self.fit_beta = None
        self.fit_done = False
        self.set_defaults()
        if init_fit:
            self.grey_body_fit()
            self.log.debug(self.info())

    def __repr__(self):
        return "Sed(source_name={}, coordinate={}, )"

    def __str__(self):
        return self.info()

    def info(self):
        string =  'Source Name: {}\n'.format(self.source_name)
        string += 'Coordinate: {}\n'.format(self.coordinate)
        string += 'flux_array: {}\n'.format(self.flux_array)
        if self.fit_done:
            string += 'Temperature Fit: {}\n'.format(self.fit_temperatures)
            string += 'Mass Fit: {}\n'.format(self.fit_masses)
            string += 'Beta Fit: {}\n'.format(self.fit_beta)
            string += 'Chisq: {}\n'.format(self.fit_chisq)
        if not self.fit_done:
            sting = 'SED not fitted yet'
        return string

    def set_defaults(self):
        # Set a default guess for the input parameter
        # needed to run a greybody Fit.
        init_guess_temperature = self.temperature_guess or [20, 50, 100, 150, 200]
        init_guess_masses = self.mass_guess or [1e5, 1e2, 1e1 ,1e-1, 1e-2]
        init_beta_guess = self.beta_guess or [2.]
        self.temperature_guess = []
        self.mass_guess = []
        for x in range(self.number_components):
            self.temperature_guess += [init_guess_temperature[x]]
            self.mass_guess += [init_guess_masses[x]]
        self.beta_guess = init_beta_guess
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

    def set_initial_guess(self, temperature_guess, mass_guess, beta_guess):
        init_guess_temperature = temperature_guess
        init_guess_masses = mass_guess
        init_beta_guess = beta_guess
        self.temperature_guess = []
        self.mass_guess = []
        if len(init_guess_temperature) < self.number_components:
            self.log.error("Not enough temperature guesses entered.")
            raise SystemExit
        if len(init_guess_masses) < self.number_components:
            self.log.error("Not enough mass guesses entered.")
            raise SystemExit
        for x in range(self.number_components):
            self.temperature_guess += [init_guess_temperature[x]]
            self.mass_guess += [init_guess_masses[x]]
        self.beta_guess = beta_guess
        self.p1 = [self.temperature_guess, self.mass_guess, self.beta_guess]
        self.p2 = None

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
                 linewidth=0.5, x_range='normal'):
        '''Plot a multi component greybody model.

        Parameters
        ----------

        nu_or_lambda :
           plot against frequency ``'nu'`` or wavelenght ``'lambda'``
        kappa :
            The kappa to use. ``'easy'`` or ``'Kruegel'``. Please refer
            to :py:func:`functions.astroFunctions.greyBody` for more
            information.
        xRange : PLEASE ADD DESCRIPTION
        linewidth : float
            The linewidth of the plotted lines. Default to 0.5.
        color : matplotlib conform color
            the color of the plotted lines. Default to ``'black'``.
        '''
        if not self.fit_done:
            raise ValueError('Could not plot the data. The data has not been '
                             'fitted yet.')
            sys.exit()
        if self.p2 ==  None:
            pass
        if x_range == 'LTIR':
        # Plot the SED in the range of the determination
            # of the L_TIR: 3-1100 micron
            xmin =  3e-6# micron
            xmax =  4000e-6 # micron
            # conversion to frequency in GHz
            xmin = const.c/xmax/1e9
            xmax = const.c/xmin/1e9
            step = 0.1

        if x_range == 'normal':
            # arbitrary range definition
            xmin = 1e-2
            xmax = 3e5
            step = 0.5
        if type(x_range) == list:
            xmin = x_range[0]
            xmax = x_range[1]
            if len(x_range) < 3:
                step = 0.1
            else:
                step = x_range[2]
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
            y=x.copy()
            modelLambda =model.copy()
            greyLambda = []
            for i in range(len(grey)):
                greyLambda += [grey[i].copy()]
            for i in range(len(x)):
                y[i]=(const.c/(x[len(x)-i-1]*1e9))/1e-6
            for i in range(len(model)):
                modelLambda[i]=model[len(model)-i-1]
                for j in range(len(greyLambda)):
                    greyLambda[j][i] = grey[j][len(grey[j])-i-1]
            x=y
            model =modelLambda
            grey = greyLambda
        plt.loglog(x, model, ls='-', color=color, label='_nolegend_', lw=0.5,
                   marker='')
        linestyles = [':','-.','-']
        j=0
        for i in grey:
            plt.loglog(x, i, color=color, ls=linestyles[j], lw=0.5, marker='')
            j+=1

    def create_figure(self, save=True, plotLegend=False,
                      color=['black'], marker=['x'], title=None, x_label=None,
                      y_label=None, nu_or_lambda='nu', fontdict=None,
                      textStringLoc=[1,1], lineWidth=0.5,
                      kappa='easy', x_range='normal', prefix="./", **kwargs):
        r""" Creates a quick preview of the loaded SED. TODO: extend
        documentation.
        """
        fig1 = plt.figure()
        fig1ax1 = fig1.add_subplot(111)
        textString = ''
        for i in range(len(self.p2[0])):
            textString += ('T=' + str('%1.1f' % self.p2[0][i]) +
                           ' K\nM=' + str("%1.2e" % self.p2[1][i]) + ' Msun\n')

        if len(self.p2[0])==2:
            textString+= 'N1/N2 = '+str('%i'%(self.p2[1][0]/self.p2[1][1]))+'\n'
        textString += ('beta = ' + str("%1.2f" % self.p2[2][0]) +
                       '\nchi$^2$ =' + str("%1.2f" % self.fit_chisq) + '\n')
        # sets the limits of the plot Page
        plotSize = 0.9 # how much in percentace should the plotpage be larger
                       # than the plotted values?
        if nu_or_lambda=='nu':
            xLimNu = [min(self.flux_array[0]) -
                      min(self.flux_array[0]) * plotSize,
                      max(self.flux_array[0]) + max(self.flux_array[0])]

        if nu_or_lambda == 'lambda':
            newSelf_Flux_Array = []
            for i in self.flux_array[0]:
                print i
                newSelf_Flux_Array += [astro_functions.frequency_to_wavelength(i)]

            self.flux_array[0] =newSelf_Flux_Array
            xLimNu = [min(self.flux_array[0]) -
                      min(self.flux_array[0]) * plotSize * 2,
                      max(self.flux_array[0]) +
                      max(self.flux_array[0]) * plotSize]
        ylim = [min(self.flux_array[1]) - min(self.flux_array[1]) * plotSize,
                max(self.flux_array[1]) +
                max(self.flux_array[1]) * plotSize / 2]
        # makes the plot page squared; TBD not really square yet
        fig1ax1.set_xlim(xLimNu[0],xLimNu[1])
        fig1ax1.set_ylim(ylim[0],ylim[1])

        # PLots the model given in self.p2
        self.plot_sed(axes=plt.gca(), nu_or_lambda=nu_or_lambda, color='black',
                      linewidth=0.5, x_range=x_range)

        markersize =7
        #Plotting the data points
        parted = [1]
        if len(parted)==1:
            fig1ax1.errorbar(self.flux_array[0], self.flux_array[1],
                         yerr=self.flux_array[2], fmt='o', marker='p',
                         mfc='None', mew=0.5, mec='#00ffff', ms=markersize,
                         color='black', lw=lineWidth)
        else:
            for i in range(len(parted)):
                if i == 0:
                    fig1ax1.errorbar(self.flux_array[0][0:parted[i]],
                                 self.flux_array[1][0:parted[i]],
                                 yerr=self.flux_array[2][0:parted[i]],
                                 fmt=marker[i],
                                 marker=marker[i], mfc='None', label=label[i],
                                 mew=0.5, mec=color[i], ms=markersize,
                                 color=color[i], lw=lineWidth)
                else:
                    fig1ax1.errorbar(self.flux_array[0][parted[i-1]:parted[i]],
                                self.flux_array[1][parted[i-1]:parted[i]],
                                yerr=self.flux_array[2][parted[i-1]:parted[i]],
                                fmt=marker[i], marker=marker[i],
                                mfc='None', label=label[i], mew=0.5,
                                mec=color[i], ms=markersize, color=color[i],
                                lw=lineWidth)

        # setting up legend,title, xlabel.
        if plotLegend == 'yes':
            fontdict={'size':'11'}
            plt.legend(loc='upper right', numpoints=1, fancybox=False,
                      prop=fontdict, markerscale=1)
        fontdict={'size':'17'}
        plt.text(90,0.2,s=textString, fontdict=fontdict, alpha=0.4)
        if title:
            fig1ax1.title(title)
        if x_label:
            fig1ax1.xlabel(x_label)
        if y_label:
            fig1ax1.ylabel(y_label)
        fig1ax1.axis([xLimNu[0],xLimNu[1],ylim[0],ylim[1]])
        if save:
            fig1.savefig(prefix + self.source_name + '_SED.eps', dpi=None,
                         facecolor='w', edgecolor='w',orientation='portrait',
                         papertype='a5', format='eps',bbox_inches='tight')
