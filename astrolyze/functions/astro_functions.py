# Copyright (C) 2012, Christof Buchbender
# BSD Licencse
import math
import sys
import os
import scipy

from scipy.optimize import leastsq as least
from numpy import asarray, mean, std, where, exp, log, sqrt, arange
from copy import deepcopy as copy
import random as rnd

from generaltools import log_tools

import astrolyze.functions.constants as const
import astrolyze.functions.units

USER = os.getenv("USER")

log = log_tools.init_logger(
    directory="/home/{}/.astrolyze/".format(USER),
    name="astrolyze"
)

def black_body(x, T, nu_or_lambda='nu'):
    r"""
    Calculation of the flux density of a black body at a temperature T and a
    wavelenght/frequency x.

    Parameters
    ----------
    x : float or numpy array
        wavelength [GHz] or frequency [micron];  specify type in nu_or_lambda
    T : float [Kelvin]
        Temperature of the black_body
    nu_or_lambda : string
        Specify whether x is a frequency :math:`\nu` ``'nu'`` or a wavelenght
        :math:`\lambda` ``'lambda'``; default is ``'nu'``.

    Returns
    -------
        Flux density in Jansky : float [Jy]

    Notes
    -----
    This functions resembles the following formulas for input in frequency:

    .. math::

        B_{\nu} = \frac{2 h \nu^3}{c^2} (e^{\frac{h \nu}{k T}} - 1)^{-1}

    and for input in wavelenght:

    .. math::

        B_{\lambda} = \frac{2 h c^2}{\lambda^5} (e^{\frac{h \lambda}{k T}} -
        1)^{-1}

    Both formulas are scaled by 1e26, thus returning the flux in Jansky.

    Examples
    --------
    The function works with linear numpy arrays. Thus the black_body can be
    evaluated at many points at the same time. Using matplotlib it can
    also be plotted:

    .. plot::
       :include-source:

       import numpy as np
       import matplotlib.pyplot as plt
       import astrolyze.functions.astro_functions as astFunc

       frequency_range = np.arange(1e4, 1e7, 1e4)
       temperature_1 = 6000
       temperature_2 = 12000  # Kelvin
       blackbody_1 = astFunc.black_body(frequency_range, temperature_1)
       blackbody_2 = astFunc.black_body(frequency_range, temperature_2)
       figure = plt.figure()
       axis = figure.add_subplot(111)
       pl = axis.loglog(frequency_range, blackbody_1, label='T = 6000 K')
       pl = axis.loglog(frequency_range, blackbody_2, label='T = 12000 K')
       pl = axis.legend()
       plt.savefig('black_body.eps')

    """
    if nu_or_lambda == 'nu':
        return ((2 * const.h * ((x * 1e9) ** 3) / (const.c ** 2)) *
                (1 / (math.e ** ((const.h * x * 1e9) / (const.k *
                float(T))) - 1)) / 1e-26)
    if nu_or_lambda == 'lambda':
        return ((2 * const.h * const.c ** 2 / ((x * 1e-6) ** 5)) *
                (1 / (math.e ** ((const.h * const.c) /
                (x * 1e-6 * const.k * float(T))) - 1)) /
                1e-26)


def grey_body(p, x, nu_or_lambda='nu', kappa='Kruegel', distance=840e3):
    r""" Calculation of the flux density in Jansky of a grey_body under
    assumption of optically thin emission. Please see Notes below for an
    detailed description assumptions and equations used.

    Parameters
    ----------
    p : list
        List of the parameters defining a grey_body, being Temperature [K],
        column density or mass (dependent on the kappa used) and the grey_body
        slope index beta, respectively (refer to notes for more information):
        p = [T, N, beta]
    x : float or numpy array
        Wavelength [GHz] or frequency [micron];
        specify type in nu_or_lambda
    kappa : string
        Chooses the dust extinction coefficient to use:
            * ``"easy"`` -> kappa = nu^beta; tau = N * kappa
            * ``"Kruegel"`` -> kappa = 0.04*(nu/250Ghz)^beta;
              tau = M/D^2 * kappa

        Please refer to Notes below, for further explanation.
    distance : float
        The distance to the source that is to be modeled if kappa
        ``"Kruegel"`` is used.

    Other Parameters
    ----------------
    nu_or_lambda : string
        Specify whether x is a frequency :math:`\nu` ``'nu'`` or a wavelength
        :math:`\lambda` ``'lambda'``; default is ``'nu'``. if lambda the input
        converted to a frequency in [GHz].

    Notes
    -----
    The general equation for a grey_body is:

    .. math::

        S(x, \tau) = (black_body(x, T)) * [1 - e^\tau] \Omega

    describing the flux coming from an solid angle
    :math:`\Omega` while :math:`\tau` is:

    .. math::

        \tau_{\nu} = \frac{ \kappa_d(\nu) * M_{dust}}{D^2 \Omega} .

    Here we assume optically thin emission and a source filling factor of
    unity. This simplifies the equation of the grey_body to:

    .. math::

        S(x, \tau) = \tau * (black_body(x, T))

    This script supports two versions of the dust extinction coefficient.::
    A simple version without a lot of physics put into, kappa = ``'easy'``
    which defaults to the following grey_body equation:

    ..  math::

        S(x, \tau) = N * x^{\beta} * black_body(x,T) ,

    with N being a column density scaling factor.

    The second version, kappa = ``'Kruegel'`` uses the dust extinction
    coefficient reported in [KS] which renders the used equation to:

    ..  math::

        \kappa = 0.04 * (\frac{x\,[GHz]}{250\,GHz})^{\beta}

        S_{\nu} = M[kg] / D^2[m^2] * \kappa * black_body(x,T) .

    Examples
    --------
    The same examples as for :func:`black_body` apply.

    References
    ----------
    .. [KS] Kruegel, E. & Siebenmorgen, R. 1994, A&A, 288, 929

    """
    T, N, beta = p
    # Switch to choose kappa
    if kappa == 'easy':
        # Here N is an arbitrary fitParameter with no direct physical meaning.
        kappa = (x * 1e9) ** beta
        tau = kappa * N
    if kappa == 'Kruegel':
        # Be careful, for kappa='Kruegel' the start Parameter N is actually
        # the Mass of the specific component in sun masses.
        kappa = 0.04 * ((x * 1e9 / 250e9) ** beta)  # [KS]
        distance_in_m = distance * const.parsec_in_m  # M33.distance is 840e3 in
                                                         # parsec. pcInM from
                                                         # constants
        D2 = distance_in_m ** 2  # calculate the Distance squared [m^2]
        M = N * const.m_sun  # Convert the Mass which is given as a start
                          # Parameter in Msun to kg
        tau = kappa * (M / D2)
    if nu_or_lambda == 'lambda':
        # Conversion of x in wavelenght [micron] to frequency [GHz].
        x = const.c / (x * 1e-6) / 1e9
        nu_or_lambda == 'nu'
    if nu_or_lambda == 'nu':
        return tau * black_body(x, T, nu_or_lambda)
#    if nu_or_lambda == 'lambda':
#        # Not used anymore
#        print 'Warning: This part of the script is not up to date'
#        return N * ((c/(x * 1e-6))**beta) * black_body(x,T,nu_or_lambda)

def multi_component_grey_body(pMulti, x, nu_or_lambda='nu', kappa='Kruegel'):
    r"""
    Combines multiple grey_body functions and returns the flux density in
    Jansky for the input frequency/wavelenght.

    Parameters
    ----------
    pMulti : nested lists
        Similar to p from
        :py:func:`grey_body` but the three entries are
        lists, i.e.::
        pMulti = [[T1, T2, T3, ...Tn], [N1, N2, N3,...Nn], [beta]]
    x : float or numpy array
        frequency [micron] (or wavelenght **Not maintained**, specify type in
        nu_or_lambda)

    Returns
    -------
    sum(snu) : float
        All dust components summed.
    snu :
        A list with the fluxes of the individual components.

    See Also
    --------
    black_body, grey_body

    Notes
    -----
    Only one common beta for all components can be used. May be expanded to
    mutliple betas if needed.

    Examples
    --------
    Same as for black_body, but all returned grey_bodies may be plotted.
    """
    T, N, beta = pMulti
    if type(beta) == list:
        beta = beta[0]
    snu = []
    for i in range(len(T)):
        pOne = [T[i], N[i], beta]
        snu += [grey_body(pOne, x, nu_or_lambda, kappa=kappa)]
    snu = asarray(snu)
    return snu.sum(axis=0), snu


def grey_body_fit(data, start_parameter, nu_or_lambda='nu', fit_beta=False,
                fix_temperature=False, rawChiSq=None, kappa='Kruegel',
                residuals=False, iterations=1e9):
    r"""
    This function fits a multi component grey body model to an observed SED for
    the optical thin case.

    Parameters
    ----------
    data : array
        The obseved data. Array of shape(3, x) first row has to be the X values
        (Frequency in [GHz]) of the measurements, second row the Y values
        (Flux [Jy]), and the third row  the Z values the errors on the
        fluxes i.e.:
        data = array([[X1, X2, X3, ...], [Y1, Y2, Y3,...], [Z1, Z2,
        Z3, ...]])
    start_parameter : array
        Array of a first guess of the parameters of the grey_body components.
        The number of components is arbitrary.
        start_parameter = [[T1, T2, T3,...], [N1, N2, N3, ...], beta]

    fit_beta : True or False
        If True Beta is allowed to vary. Default is False.

    fix_temperature : True or False
        If True the Temperature is fixed allowed to vary.

    rawChiSq :
        if None the function gives the reduced chisq Value. If True the
        function gives chisq without dividing it by the dof

    Returns
    -------
    p2 : list
        The final grey_body parameters that reduce the least squares for the
        given dataset.
    chisq/rawChiSq :
        chisq is reduced chisq with degrees of freedom:
        dof= #dataPoints-#freeFitParameters-1

    Other Parameters
    ----------------
    nu_or_lambda : string
        Specify whether x is a frequency :math:`\nu` ``'nu'`` or
        a wavelenght :math:`\lambda` ``'lambda'``; default is ``'nu'``.::
        **Don't** use ``'lambda'`` as this part of the
        :py:func:`grey_body` is not up-to-date.

    See Also
    --------
    scipy.optimize.leastsq: This function is used to perform the least squares
    fit.
    multi_component_grey_body, grey_body, black_body: Defining the function to
    be fittet to the data.

    Notes
    -----
    A one component fit has four free parameters if beta is allowed to vary or
    three if beta is fixed (one more than parameters to fit). Each additional
    component adds two more free paramters to fit.
    Assure that:
    number of data points > number of free parameters.
    """
    # Distinguish between the different options and build the parameter list
    # needed by optimize.leastsq()
    if not fit_beta:
        if not fix_temperature:
            p = start_parameter[0] + start_parameter[1]
            beta = start_parameter[2]
        if fix_temperature:
            p = start_parameter[1]
            Temps = start_parameter[0]
            beta = start_parameter[2]
    if fit_beta:
        if not fix_temperature:
            p = start_parameter[0] + start_parameter[1] + start_parameter[2]
        if fix_temperature:
            p = start_parameter[1] + start_parameter[2]
            Temps = start_parameter[0]

    def _err(p, x, y, y_error, nu_or_lambda):
        """
        The function to be minimized by scipy.optimize.leastsq. It returns the
        difference between the measured fluxes, y, and the calculated fluxes
        for the parameters p of the SED at the given frequency x.

        Parameters
        ----------
        p : list
            same start start_parameter as in grey_body_fit
        x and y :
         The two rows, data[0] and data[1] of the data variable from
         grey_body_fit.
        y_error : The absolute error on the flux measurement.

        Other Parameters
        ----------------
        nu_or_lambda : string
            Specify whether x is a frequency :math:`\nu` ``'nu'`` or
            a wavelenght :math:`\lambda` ``'lambda'``; default is ``'nu'``.::
            **Don't** use ``'lambda'`` as this part of the
            :py:func:`grey_body` is not up-to-date.

        Notes
        -----
        This function caluclates the residuals that are minimized by
        :func:`leastsq`, thus it calculates:

        .. math::

            (model-measurement)/\sigma

        Note the :math:`\chi^2` is defined as:

        .. math::

            \chi^2 = \frac{1}{N} \sum{(model-measurement)/\sigma}

        :warning:`Formula has to be checked.`
        """
        if not fit_beta:
            if not fix_temperature:
                numberOfComponents = (len(p)) / 2
                pMulti = [list(p[0:numberOfComponents]),
                          list(p[numberOfComponents:numberOfComponents * 2])]
                pMulti += [beta]
            if fix_temperature:
                pMulti = [Temps, p, beta]

        if fit_beta:
            if not fix_temperature:
                numberOfComponents = (len(p) - 1) / 2
                pMulti = [list(p[0:numberOfComponents]),
                          list(p[numberOfComponents:numberOfComponents * 2]),
                          p[len(p) - 1]]

            if fix_temperature:
                numberOfComponents = len(p) - 1
                pMulti = [Temps, list(p[0:numberOfComponents]), p[len(p) - 1]]

        return (((multi_component_grey_body(pMulti, x, nu_or_lambda,
                kappa=kappa)[0]) - y) / y_error)
    # The actual fit
    # maxfev : Number of integrations
    # full_output: Additional informations
    # args: X and Y data
    p2, cov, info, mesg, success = least(_err, p, args=(data[0], data[1],
                                         data[2], nu_or_lambda),
                                         maxfev=int(iterations), full_output=1)
    # return of the optimized parameters
    dof = len(data[0]) - len(p) - 1  # degrees of Freedom
    chisq = sum(info['fvec'] * info['fvec']) / dof  # see above
    rawchisq = sum(info['fvec'] * info['fvec'])

    if not fit_beta:
        if not fix_temperature:
            numberOfComponents = (len(p)) / 2
            p2 = [list(p2[0:numberOfComponents]),
                  list(p2[numberOfComponents:numberOfComponents * 2]), beta]
        if fix_temperature:
            numberOfComponents = len(p)
            p2 = [Temps, p2, beta]

    if fit_beta:
        if not fix_temperature:
            numberOfComponents = (len(p) - 1) / 2
            p2 = [list(p2[0:numberOfComponents]),
                  list(p2[numberOfComponents:numberOfComponents * 2]),
                  [p2[len(p) - 1]]]
        if fix_temperature:
            numberOfComponents = (len(p) - 1)
            p2 = [Temps, list(p2[0:numberOfComponents]), [p2[len(p) - 1]]]
    if residuals:
        return p2, chisq, info['fvec']
    if rawChiSq == None:
        return p2, chisq
    if rawChiSq:
        return p2, rawchisq


def _grid_fit(data, beta, nu_or_lambda='nu', fit_beta=False, rawChiSq=False,
             kappa='Kruegel'):
    r"""
    Different approach to find the best fitting SED using certain ranges of
    temperatures and masses to avoid unreasonably high values.

    Not sure about the functionality and the code is badly written with the
    temperature and mass ranges hard coded. Thats why its not public.
    Once approved it may be made public again.

    Parameters
    ----------
    data : array
        Same format as in grey_body_fit.
    beta :
        The beta value. If fit_beta=True this is varied in the fits.
    kappa :
        As in :func:`greybody`.

    Other Parameters
    ----------------
    fit_beta: True or False
        Steers if beta is fittet or not.
    rawChisq: True or False
        Returns chi square without normalisation, or not.
    """
    def _err(p, x, y, y_error, nu_or_lambda):
        """
        The function to be minimized by scipy.optimize.leastsq. It returns the
        difference between the measured fluxes, y, and the calculated fluxes
        for the parameters p of the SED at the given frequency x.

        Parameters
        ----------
        p : list
            same start start_parameter as in grey_body_fit
        x and y :
         The two rows, data[0] and data[1] of the data variable from
         grey_body_fit.
        y_error : The absolute error on the flux measurement.

        Other Parameters
        ----------------
        nu_or_lambda : string
            Specify whether x is a frequency :math:`\nu` ``'nu'`` or
            a wavelenght :math:`\lambda` ``'lambda'``; default is ``'nu'``.::
            **Don't** use ``'lambda'`` as this part of the
            :py:func:`grey_body` is not up-to-date.

        Notes
        -----
        This function caluclates the residuals that are minimized by
        :func:`leastsq`, thus it calculates:

        .. math::

            (model-measurement)/\sigma

        Note the :math:`\chi^2` is defined as:

        .. math::

            \chi^2 = \frac{1}{N} \sum{(model-measurement)/\sigma}

        :warning:`Formula has to be checked.`

        This functions is different from the _err functions in grey_body_fit
        because the start parameters are handles differently.
        """
        if not fit_beta:
            pMulti = [T, p, beta]
            print pMulti
        if fit_beta:
            pMulti = [Temps, p]
            print pMulti
        return (((multi_component_grey_body(pMulti, x, nu_or_lambda, kappa)[0])
                - y) / y_error)
    beta = [beta]
    T1 = arange(5, 70, 1)
    T2 = arange(20, 100, 1)
    chisqList = {}
    chisqList1 = []
    for i in T1:
        for j in T2:
            T = [i, j]
            N = [1e2, 1]
            if not fit_beta:
                p = N
            if fit_beta:
                p = N + beta
            p2, cov, info, mesg, success = least(_err, p, args=(data[0],
                                                 data[1], data[2],
                                                 nu_or_lambda),
                                                 maxfev=int(1e9),
                                                 full_output=1)
            dof = len(data[0]) - len(p) - 1  # Degrees of Freedom
            chisq = sum(info['fvec'] * info['fvec']) / dof  # see above
            rawchisq = sum(info['fvec'] * info['fvec'])
            chisqList[chisq] = [[i, j], p2]
            chisqList1 += [chisq]
    chisq = sorted(chisqList)[0]
    T, p2 = chisqList[sorted(chisqList)[0]]
    if not fit_beta:
        numberOfComponents = len(p)
        p2 = [T, p2, beta]
    if fit_beta:
        numberOfComponents = (len(p) - 1)
        p2 = [Temps, list(p2[0:numberOfComponents]), [p2[len(p) - 1]]]
    if not rawChiSq:
        return p2, chisq
    if rawChiSq:
        return p2, rawchisq
    return sorted(chisq)[0]


def LTIR(p2, kappa='Kruegel', xmin=3., xmax=1100., beamConv=True,
         distance=847e3, unit='JyB'):
    r"""
    Integration of a multi-component greybody model.

    Parameters
    ----------
    p2 : list
        The parameters defining the multi-component greybody model. Same format
        as p in
        :py:func:`astrolyze.functions.astroFunctions.multi_component_grey_body`
    kappa : string
        The dust extinction coefficient used to describe the greybodies. See:
        py:func:`grey_body`
    xmin, xmax : float
        The integration range in units of micron. Defaults to 3 -- 110 micron.
        The definition of LTIR from [DA]
    beamConv : True or False
        For units in Lsun the code is not well written. Hardcoded conversion
        between an 28" and 40" beam. !! CHANGE !!
    unit : string
        If ``'Lsun'`` the returned integrated flux is in units of solar
        luminosities (erg s^-1). For this a distance is needed. If ``'JyB'``
        the units are Jy/beam; distance is not used.

    Notes
    -----
    Needs some work to be generally usable. For units in Jy/beam the code seems
    to be safe.

    References
    ----------
    .. [DA] Dale et al. 2001; ApJ; 549:215-227
    """
    #Convert the input xmin/xmax from micron to m and to frequency in
    #GHz
    xmin = const.c / (xmin * 1e-6) / 1e9        # GHz
    xmax = const.c / (xmax * 1e-6) / 1e9
    step = 0.1
    x = arange(floor(xmax), floor(xmin), step)
    # multi_component_grey_body needs input (x) in GHz
    model, grey = multi_component_grey_body(p2, x, 'nu', kappa)
    # integrate the SED
    LTIR1 = sum(model) * (step * 1e9)  # Jy/Beam*Hz
    if unit == 'Lsun':
        # conversion factor between 40 and 28 arcsc beam
        beamsize_conversion = (28. ** 2) / (40 ** 2)
        if beamConv:
            conv = (units.JanskyToErgs_m * units.Int2Lum(distance,
                    cm_or_m='m') * beamsize_conversion)
        # convert to erg s-1 /beam in terms of 1e6*Lsun
        conv = Jy_to_FluxDensityInErg_s * units.Int2Lum(distance, cm_or_m='m')
        LTIR1 = LTIR1 * conv / luminositySun / 1e6
    if unit == 'JyB':
        pass
    return LTIR1


def generate_monte_carlo_data_sed(data):
    """
    MonteCarlo Simulation of a set of flux measurements, assuming that the
    measurement data follows a gauss distribution.

    This function makes use of the :func:`random.gauss` function to generate a
    data point from a gauss distribution, that has a mean equal to the Flux
    measurement and a standard deviation correponding to the error of the
    measurement.

    Parameters
    ----------
    data : array
         Same format as in grey_body_fit function:
            data= [[x1, x2, x3, ...][y1, y2, y3, ...][z1, z2, z3, ...]]

         with x = wavelenght/frequency, y = flux, z = error on flux.

    Returns
    -------
    newData : array in same format as data.
        The Monte-Carlo simulated measurement.

    See Also
    --------
    random.gauss

    """
    newData = copy(data)
    for i in range(len(data[1])):
        newData[1][i] = rnd.gauss(data[1][i], data[2][i])
    return newData


def grey_body_monte_carlo(p, data, iterations):
    """
    Function to evaluate the errors in the parameters fitted with the
    grey_body_fit function.

    It uses Monte Carlo Simulated data (from
    :func:`generate_monte_carlo_data_sed`) and performs a fit to this new data
    giving back the results of the fit parameters.

    Parameters
    ----------
    p : list
        The parameters defining the multi component grey_body model to be
        fitted. Same format as p in :py:func:`multi_component_grey_body`
    data : array
        The actual measured data of the SED, same format as for
        :py:func:`grey_body_fitFunction`
    iterations : int
        Number of times new data is generated and fitted.

    Returns
    -------
    string :
        Containing the mean, standard deviation of the fit parameters, ready
        to print out.
    betaTlist : List of all fit results. Name misleading since it may not
        include the beta.
    """
    # Define the variables that store the MOnte Carlo T and N values
    TList = []
    NList = []
    chisqList = []
    betaList = []
    for i in range(len(p[0])):
        TList += [[]]
        NList += [[]]
    if len(p[0]) == 1:
        betaTList = [[], [], [], [], []]
    if len(p[0]) == 2:
        betaTList = [[], [], [], [], [], []]
    string = ''

    for i in range(0, iterations):
        print i + 1, '\\', iterations
        sys.stdout.flush()
        MCData = generate_monte_carlo_data_sed(data)
        if len(p[0]) == 1:
            p2, chisq = grey_body_fit(MCData, p, fit_beta=True,
                                   fix_temperature=False)
        else:
            p2, chisq = grey_body_fit(MCData, p, fit_beta=False,
                                   fix_temperature=False)
        chisqList += [chisq]
        x = 0
        if len(p[0]) == 1:
            betaTList[3] += [MCData]
            betaTList[4] += [p2]
        if len(p[0]) == 2:
            betaTList[4] += [MCData]
            betaTList[5] += [p2]
        for i in range(len(p[0])):
            if len(p[0]) == 1:
                betaTList[0] += [p2[0][i]]
                betaTList[1] += [p2[2][i]]
                betaTList[2] += [p2[1][i]]
            if len(p[0]) == 2:
                betaTList[i] += [p2[0][i]]
                betaTList[i + 2] += [p2[1][i]]
                #if float(p2[0][0]) < 100:
                    #if float(p2[0][1]) < 100:
                TList[i] += [p2[0][i]]
                NList[i] += [p2[1][i]]
                if x == 0:
                    betaList += p2[2]
            else:
                #if float(p2[0][i]) < 100:
                TList[i] += [p2[0][i]]
                NList[i] += [p2[1][i]]
                if x == 0:
                    betaList += p2[2]
            x = x + 1
    betaList = asarray(betaList)
    for i in range(len(p[0])):
        string += ('T' + str(i + 1) + ': ' +
                   str("%1.2f" % mean(asarray(TList[i]))) +
                   ' +/- ' + str("%1.2f" % std(asarray(TList[i]))) + '\n')
        string += ('N' + str(i + 1) + ': ' +
                   str("%1.2e" % mean(asarray(NList[i]))) +
                   ' +/- ' + str("%1.2e" % std(asarray(NList[i]))) + '\n')
    string += ('Beta: ' + str("%1.2f" % mean(betaList)) + ' +/- ' +
               str("%1.2f" % std(betaList)) + '\n')
    string += 'Number of Fits ' + str(len(TList[0])) + '\n'
    if len(p[0]) == 1:
        return string, betaTList
    else:
        return string, betaTList


def line(p, x):
    r"""
    Line `y = m*x + b` equation. Returns y value at point x.

    Parameters
    ----------
    p : list
        Contains the slope and the y-axis intersection of the line [m, b].

    Returns
    -------
    y : value of y corresponding to x.
    """
    return p[0] * x + p[1]


def anti_line(p, y):
    r"""
    Inverse of a line returning the x value corresponding to a y value, i.e.
    `x = y/m - b`.

    Parameters
    ----------
    p : list
        Contains the slope and the y-axis intersection of the line [m, b].

    Returns
    -------
    y : value of x corresponding to y.
    """
    return y / p[0] - p[1]


def linear_error_function(p, x, y, y_error, x_error):
    """
    Error function, i.e. residual from the measured value, which has to be
    minimised in the least square fit taking X and Y Error into account.

    Parameters
    ----------

    p : list
        Same as in :func:`line` and :func:`anti_line`.
    x : float or list
        x measurements. Data.
    y : float or list
        y measurements. Data.
    x_error : float or list
        The x measurment errors.
    y_error : float or list
        The y measurment errors.
    """
    if x_error.all():
        return sqrt(((line(p, x) - y) / y_error) ** 2 + ((anti_line(p, y) - x)
                    / x_error) ** 2)
    if not x_error.all():
        return sqrt(((line(p, x) - y) / y_error) ** 2)

def line_fit(p, x, y, y_error, x_error=False, iterations=10000):
    """
    Linear Fit to data, taking either errors in y or both in x and y into
    account.

    Parameters
    ----------

    p : list
        Containg slope (m) and y-axis intersection (b) p=[m, b]. Same as in
        :func:`line` and :func:`antiline`.
    x : float or list
        x measurements. Data.
    y : float or list
        y measurements. Data.
    y_error : float or list
        The y measurment errors.
    x_error : float or list
        The x measurment errors. If unset only errors in y are taken into
        account.
    """
    p1, cov, info, mesg, success = least(linear_error_function, p,
                                         args=(x, y, y_error, x_error),
                                         maxfev=int(iterations), full_output=1)
    dof = len(x) - len(p) - 1  # degrees of Freedom
    chisq = sum(info['fvec'] * info['fvec']) / dof
    Error = std(info['fvec'])
    return p1, chisq, Error


def analytic_linear_fit(x, y, x_error, y_error):
    r"""
    This function resembles the analytical solution following chaper 8 from
    [TA].

    Parameters
    ----------

    x : float or list
        x measurements. Data.
    y : float or list
        y measurements. Data.
    y_error : float or list
        The y measurment errors.
    x_error : float or list
        The x measurment errors. If unset only errors in y are taken into
        account.

    Notes
    -----

    Without errors the following holds:

    .. math::
        y = A + B x

        A = \frac{\Sigma(x^2) \cdot \Sigma(y) - \Sigma(x) \cdot
        \Sigma(x \cdot y)}{\Delta}

        B = N \frac{\Sigma(x \cdot y) - \Sigma (x) \cdot \Sigma(y)}{\Delta}

        \Delta = N \cdot \Sigma(x^2) - (\Sigma(x))^2

    .. warning:: This has to be checked.

    References
    ----------
    .. [TA] "An introduction to the study of uncertainties in physical
        measurement" by John R. Taylor.
    """
    # first calculate a least squares fit ignoring the errors since B is
    # needed for the more complex issue including errors
    sumX = sum(x)
    sumY = sum(y)
    sumXY = sum(x * y)
    sumXSq = sum(x ** 2)
    N = len(x)
    Delta = N * sumXSq - (sumX) ** 2
    A = (sumXSq * sumY - sumX * sumXY) / Delta
    B = (N * sumXY - sumX * sumY) / Delta
    print 'm = ' + '%1.2f' % A + ' b = ' + '%1.2f' % B
    # now make use of the idea of a equivalent error only in y defined by
    # equivalentError =  sqrt(y_error**2+(B*x_error)**2)
    equivalentError = y_error  # sqrt(y_error**2+(B*x_error)**2)
    # and use wheighted least square fit see definition for A and B and
    # their errors below
    weight = 1. / (equivalentError) ** 2
    sumWeightX = sum(weight * x)
    sumWeightY = sum(weight * y)
    sumWeightXY = sum(weight * x * y)
    sumWeightXSq = sum(weight * x ** 2)
    sumWeight = sum(weight)
    WeightDelta = (sumWeight * sumWeightXSq) - ((sumWeightX) ** 2)
    A = (((sumWeightXSq * sumWeightY) - (sumWeightX * sumWeightXY)) /
         WeightDelta)
    B = ((sumWeight * sumWeightXY) - (sumWeightX * sumWeightY)) / WeightDelta
    SigmaA = sqrt(sumWeightXSq / WeightDelta)
    SigmaB = sqrt(sumWeight / WeightDelta)
    Chisq = sum((y - A - B * x) ** 2 / equivalentError ** 2) / (N - 2)
    print ('m = ' + '%1.2f' % A + '+/-' + '%1.2f' % SigmaA + ' b = ' + '%1.2f'
           % B + '+/-' + '%1.2f' % SigmaB + ' chisq:' + '%1.2f' % Chisq)
    return A, SigmaA, B, SigmaB, Chisq


def generate_monte_carlo_data_line(data, errors):
    """
    This function makes a Monte Carlo Simulation of a data Set of measurements
    it uses the random.gauss() function to generate a data point
    from a gauss distribution, that has a mean equal to the measurement
    and its standard deviation corresonding to the error of the measurement.

    Parameters
    ----------

    data : list
        A list of original measurements.
    errors : list
        A list of the corresponding errors.

    Returns
    -------

    newData : array in same format as data.
        The monte carlo simulated measurement.

    See Also
    --------

    random.gauss
    """

    newData = copy(data)
    for i in range(len(data)):
        newData[i] = rnd.gauss(data[i], errors[i])
    return newData


def line_monte_carlo(p, x, y, x_error, y_error, iterations,
                              fitIterations=1e9):
    """
    Gererate an estimate of the errors of the fitted parameters determined by
    the :py:func:`line_fit` function.

    Parameters
    ----------

    p : list
        Containg slope (m) and y-axis intersection (b) p=[m, b]. Same as in
        :func:`line` and :func:`antiline`.
    x : float or list
        x measurements. Data.
    y : float or list
        y measurements. Data.
    y_error : float or list
        The y measurment errors.
    x_error : float or list
        The x measurment errors. If unset only errors in y are taken into
        account.

    Returns
    -------

    string : A string containing the results.
    BList : A list containing the fittet y-Axis intersections.
    MList : A list containing the fitted slopes.
    chisqList : A list with the chisq values.
    resultArray : Array with the mean and the standard deviations of
        slopes and y-axis intersections, i.e. [mean(M), std(M), mean(B),
        std(B)]

    See Also
    --------

    grey_body_fit, generate_monte_carlo_data_line
    """
    # Define the variables that store the MOnte Carlo B and M values y= mx+b
    BList = []
    MList = []
    chisqList = []
    FitList = []
    string = ''
    for i in range(0, iterations):
        sys.stdout.write(str(i + 1) + '\\' + str(iterations) + '\r')
        sys.stdout.flush()
        xMCData = generate_monte_carlo_data_line(x, x_error)
        yMCData = generate_monte_carlo_data_line(y, y_error)
        p2, chisq, Error = line_fit(p, xMCData, yMCData, x_error, y_error,
                                 XY_or_Y='XY', iterations=10000)
        chisqList += [chisq]
        BList += [p2[1]]
        MList += [p2[0]]
    string += ('B: ' + str("%1.2f" % mean(asarray(BList))) + ' +/- ' +
               str("%1.2f" % std(asarray(BList))) + '\n')
    string += ('M: ' + str("%1.2f" % mean(asarray(MList))) + ' +/- ' +
               str("%1.2f" % std(asarray(MList))) + '\n')
    string += ('Chi: ' + str("%1.2f" % mean(asarray(chisqList))) + ' +/- ' +
               str("%1.2f" % std(asarray(chisqList))) + '\n')
    string += 'Number of Fits ' + str(len(BList)) + '\n'
    resultArray = [mean(asarray(MList)), std(asarray(MList)),
                   mean(asarray(BList)), std(asarray(BList))]
    return string, BList, MList, chisqList, resultArray


def gauss1D(x, fwhm, offset=0, amplitude=1):
    r"""
    Calulcates 1D Gaussian.

    Parameters
    ----------
    x : float or numpy.ndarray
        the x-axis value/values where the Gaussian is to be caluclated.
    fwhm : float
        The width of the Gaussian.
    offset :
        The offset in x direction from 0. Default is 0.
    amplitude :
        The height of the Gaussian. Default is 1.

    Returns
    -------
    gauss : float or np.ndarray
        The y value for the specified Gaussian distribution evaluated at x.

    Notes
    -----
    The function used to describe the Gaussian is:

    .. math::

        f = \frac{1}{fwhm * sqrt(2 * \pi)} * e^{-1/2 (\frac{x-x0}{fwhm})^2}
    """
    gauss = 1 / 2 * ((x - offset) / fwhm) ** 2
    gauss = exp(-1 * gauss)
    gauss = 1 / (fwhm * math.sqrt(2 * math.pi)) * gauss
    return gauss


def gauss2D(x, y, major, minor, pa=0, xOffset=0, yOffset=0, amplitude=1):
    r"""
    Calculates a 2D Gaussian at position x y.

    Parameters
    ----------
    x : float or numpy.ndarray
        the x-axis value/values where the Gaussian is to be caluclated.
    y : float or numpy.ndarray
        the y-axis value/values where the Gaussian is to be calculated.

    major, minor : float
        The fwhm of the Gaussian in x and y direction.
    pa : float
        The position angle of the Gaussian in degrees. Default is 0.
    xOffset, yOffset:
        The offset in x and y direction from 0. Default is 0.
    amplitude :
        The height of the Gaussian. Default is 1.

    Returns
    -------
    gauss : float or np.ndarray
        The y value for the specified Gaussian distribution evaluated at x.

    Notes
    -----
    The function used to describe the Gaussian is :

    .. math::

        f = (amplitude * exp (-1 (a*(x-xOffset)^2 + 2*b*(x-xOffset)*(y-yOffset)
            + c*(y-yOffset)^2)))

    where:

    .. math::

        a = cos(pa)**2/(2*major**2) + sin(pa)**2/(2*minor**2) \\
        b = (-1*sin(2*pa)/(4*major**2))+(sin(2*pa)/(4*minor**2)) \\
        c = sin(pa)**2/(2*major**2) + cos(pa)**2/(2*minor**2) \\
    """
    pa = pa * math.pi / 180
    a = cos(pa) ** 2 / (2 * major ** 2) + sin(pa) ** 2 / (2 * minor ** 2)
    b = ((-1 * sin(2 * pa) / (4 * major ** 2)) + (sin(2 * pa) / (4 * minor **
         2)))
    c = sin(pa) ** 2 / (2 * major ** 2) + cos(pa) ** 2 / (2 * minor ** 2)
    gauss = a * (x - xOffset) ** 2
    gauss += 2 * b * (x - xOffset) * (y - yOffset)
    gauss += c * (y - yOffset) ** 2
    gauss = exp(-1 * gauss)
    gauss = amplitude * gauss


def degrees_to_equatorial(degrees):
    r"""
    Converts RA, DEC coordinates in degrees to equatorial notation.

    Parameters
    ----------
    degrees : list
        The coordinates in degrees in the format of: [23.4825, 30.717222]

    Returns
    -------
    equatorial : list
        The coordinates in equatorial notation, e.g.
        corresponding ['1:33:55.80', '+30:43:2.00'].
    """
    coordinate = []
    coordinate += [str(int(degrees[0] / 15)) + ':' + str(int(((degrees[0] / 15)
              - int(degrees[0] / 15)) * 60)) + ':' + "%1.2f" %
              (float(str((((degrees[0] / 15 - int(degrees[0] / 15)) * 60) -
              int((degrees[0] / 15 - int(degrees[0] / 15)) * 60)) * 60)))]
    coordinate += [(str(int(degrees[1])) + ':' +
              str(int(math.fabs(int((float(degrees[1]) - int(degrees[1])) *
              60)))) + ':' + "%1.2f" %
              (math.fabs(float(str(float(((float(degrees[1]) - int(degrees[1]))
              * 60) - int((float(degrees[1]) - int(degrees[1])) * 60)) *
              60)))))]
    return coordinate


def equatorial_to_degrees(equatorial):
    r"""
    Converts RA, DEC coordinates in equatorial notation to degrees.

    Parameters
    ----------
    equatorial : list
        The coordinates in degress in equatorial notation, e.g.
        ['1:33:55.80', '+30:43:2.00']

    Returns
    -------
    degrees : list
        The coordinates in degreees, e.g. [23.4825, 30.717222].

    Raises
    ------
    SystemExit
        If ``equatorial`` is not a list of strings in the above format.
    """
    try:
        CoordsplitRA = equatorial[0].split(':')
        CoordsplitDec = equatorial[1].split(':')
    except AttributeError:
        log.debug("equatorial_to_degrees needs a pair of RA DEC "\
                       "coordinated in equatiorial notation as input")
        raise
    if float(CoordsplitDec[0]) > 0:
        degrees = [(float(CoordsplitRA[0]) * (360. / 24) +
                  float(CoordsplitRA[1]) * (360. / 24 / 60) +
                  float(CoordsplitRA[2]) * (360. / 24 / 60 / 60)),
                  (float(CoordsplitDec[0]) + float(CoordsplitDec[1]) * (1. /
                  60) + float(CoordsplitDec[2]) * 1. / 60 / 60)]
    if float(CoordsplitDec[0]) < 0:
        degrees = [(float(CoordsplitRA[0]) * (360. / 24) +
                  float(CoordsplitRA[1]) * (360. / 24 / 60) +
                  float(CoordsplitRA[2]) * (360. / 24 / 60 / 60)),
                  (float(CoordsplitDec[0]) - float(CoordsplitDec[1]) * (1. /
                  60) - float(CoordsplitDec[2]) * 1. / 60 / 60)]
    return degrees


def calc_offset(central_coordinate, offset_coordinate, angle = 0,
                output_unit='farcsec'):
    r"""
    Calculates the offset between two coordinates.

    Parameters
    ----------
    central_coordinate : list
        The reference coordinate in degrees or equatorial.
    offset_coordinate : list
        The second coordinate, the offset will be with rescpect to
        central_coordinate.
    angle : float
        The angle in degrees, allowing rotated systems.

    Returns
    -------
    rotated_offset : list
        The offsets, rotated only if angle given.

    Notes
    -----
    This functions includes a correction of the RA offset with declination:

    .. math:

        ra_corrected = ra cos(dec)
    """
    possible_units = ['DEGREE', 'DEGREES', 'ARCMINUTE', 'ARCMINUTES', 'ARCSEC',
                      'ARCSECS']
    if output_unit.upper() not in possible_units:
        raise ValueError('Unit has to be one of the following. "' +
                         '" "'.join(possible_units).lower() + '"')
    angle =  math.radians(angle)
    central_in_degrees = equatorial_to_degrees(central_coordinate)
    offset_in_degrees = equatorial_to_degrees(offset_coordinate)
    offset = [offset_in_degrees[0] - central_in_degrees[0] ,
              offset_in_degrees[1] - central_in_degrees[1]]
    # correction for declination
    offset = [offset[0] * math.cos(math.radians(offset_in_degrees[1])),
              offset[1]]
    # Rotate the offsets.
    rotated_offset = rotation_2d(offset, angle)
    rotated_offset = asarray(rotated_offset)
    if output_unit.upper() in ['DEGREE', 'DEGREES']:
        pass
    if output_unit.upper() in ['ARCMINUTE', 'ARCMINUTES']:
        rotated_offset = rotated_offset * 60
    if output_unit.upper() in ['ARCSEC', 'ARCSECS']:
        rotated_offset = rotated_offset * 60 * 60
    return rotated_offset


def rotation_2d(coordinate, angle):
    r"""
    Implementation of the rotation matrix in two dimensions.

    Parameters
    ----------

    coordinates : list of floats
        Coordinates in the unrotated system [x, y].
    angle : float
        The rotation angle

    Returns
    -------

    [x_rotated, y_rotated]: list of floats
        Coordinates in the rotated system.
    """
    x, y = coordinate
    x_rotated = math.cos(angle) * x - math.sin(angle) * y
    y_rotated = math.sin(angle) * x + math.cos(angle) * y
    return [x_rotated, y_rotated]


def vel_to_freq_resolution(center_frequency, velocity_resolution):
    r""" Converts a velocity resolution to frequency resolution for a given
    center frequency.

    Parameters
    ----------
    center_frequency : float
        Center frequency in GHz.
    velocity_resolution :
        Velocity resolution in km/s.

    Returns
    -------
    frequency_resolution : float
        The corresponding frequency resolution in Mhz

    Notes
    -----
    Approved!

    """
    # Conversion from km/s to m/s
    velocity_resolution = velocity_resolution * 1e3
    # Conversion from GHz to Hz
    center_frequency = center_frequency * 1e9
    # Calculation of the frequency_resolution in Hz
    frequency_resolution =  (-1 * velocity_resolution *
                              center_frequency / const.c)
    # Conversion to MHz
    frequency_resolution = frequency_resolution / 1e6
    return frequency_resolution

def freq_to_vel_resolution(center_frequency, frequency_resolution):
    r""" Function to convert a frequency resolution to a velocity resolution
    for a given center frequency.

    Parameters
    ----------
    center_frequency : float
        Center frequency in GHz.
    frequency_resolution : float
        The frequency resolution in MHz.

    Returns
    -------
    velocity_resolution in km/s.

    Notes
    -----
    Uses the formula TODO v_LSR = c((nu0-nuObs)/nu0)

    Approved!
    """
    center_frequency = center_frequency * 1e9
    frequency_resolution = frequency_resolution * 1e6
    observation_frequency = center_frequency + frequency_resolution
    velocity_resolution = v_lsr(center_frequency,
                                          observation_frequency)
    # Difference between nu0 and nuObs is the velocity resolution
    return velocity_resolution

def v_lsr(center_frequency, observation_frequency):
    r""" Calculates the velocity that corresponds to a certain frequency shift
    between two frequencies.

    Parameters
    ----------
    center_frequency : float
        center_frequency in GHz
    observation_frequency : float
        The observation frequency in GHz.

    Returns
    -------
    v_lsr : float
        The velocity corresponding to the frequency shift in km/s

    Notes
    -----
    Approved!
    """
    center_frequency = center_frequency * 1e9
    observation_frequency = observation_frequency * 1e9
    v_lsr = (const.c * ((center_frequency - observation_frequency) /
                        center_frequency) / 1e3)
    return v_lsr

def redshifted_frequency(rest_frequency, v_lsr):
    r""" Calculates the sky frequency corresponding to a rest frequency for a
    source with a velocity v_lsr.

    Parameters
    ----------
    rest_frequency : float
        The frequency of the line at rest in Ghz (More often state the obvious
        :)).
    v_lsr : float
        The velocity of the source in km/s.

    Returns
    -------
    redshifted_frequency : float
       The sky frequency in GHz.

    Notes
    -----
    The formula used is:

    .. math::

        \nu_{sky} = \nu_{rest} * \frac{-1 v_{lsr}}{c + 1}

    Approved!
    """
    # Convert frequency to Hz,
    rest_frequency =  rest_frequency * 1e9
    # Convert velocity to m/s,
    v_lsr = v_lsr * 1e3
    # Calculate the sky frequency,
    redshifted_frequency = (-1. * v_lsr / const.c + 1) * rest_frequency
    # Convert to GHz.
    redshifted_frequency = redshifted_frequency / 1e9
    return redshifted_frequency


def frequency_to_wavelength(frequency):
    r"""
    Converting frequency to wavelength.

    Parameters
    ----------

    frequency : float [GHZ]

    Returns
    -------

    wavelength : float [micron]
    """
    # Conversion from GHz to Hz
    frequency = frequency * 1e9
    wavelength = const.c / frequency
    # Conversion from m to micron (mum).
    wavelength = wavelength / 1e-6
    return wavelength


def _equatorial2DegFile(inputFile):
    '''
    Old Functions if needed can be made public...
    converts equatorial coordinates to degrees
    format of input File must be:
    sourcName Ra Dec
    with Space/tab between the entries
    '''
    filein = open('positions.txt').readlines()
    coords = []
    for i in filein:
        i = i.split()
        coords+=[[i[1], i[2]]]
    print coords
    for i in coords:
        print (filein[x].split()[0], equatorial2Deg(i)[0], ',',
               equatorial2Deg(i)[1])
        x+=1


if __name__ == "__main__":
    import doctest
    doctest.testmod()
