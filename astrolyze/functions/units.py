# Copyright (C) 2012, Christof Buchbender
# BSD Licencse
import math

import astrolyze.functions.constants as const

# Constant conversion factors.
#==============> Approved !!! <==========================
WattToErgs    = 1e7  # 1W = 1e7 erg/s
ErgsToWatt    = 1e-7  # 1W = 1e-7 erg/s
JanskyToWatt  = 1e-26  # 1Jy = 1e-26 W/m2/Hz
WattToJansky  = 1e26  # 1W  = 1 Jy * m2 * Hz
ErgsToJansky_cm  = 1e23  # 1 erg/s =  1e23 Jy * cm2 * Hz * s
JanskyToErgs_cm  = 1e-23  # 1 Jy = 1e-23 erg/s/cm2/Hz
ErgsToJansky_m  = 1e19  # 1 erg/s = 1e-19 Jy * m2 * Hz * s
JanskyToErgs_m  = 1e-19  # 1 Jy = 1e-19 erg/s/m2/Hz
#==============> Approved !!! <==========================


def WmToKkms(x, resolution=0, sterad=False, ToKKms=False, m2_or_cm2='m',
             nu_or_lambda='nu'):
    '''
    Conversion between W/m2 and K km/s.

    Parameters
    ----------
    x : float
        wavelenght/frequency [GHZ].
    resolution : float
    ToKKms : True or False
        Direction of the conversion.
    sterad : True or False
        If False convert from per beam to per sterad.
    m2_or_cm2 : string
        Choose if conversion to/from W m-2 oder W cm-2. ``'m2'`` or ``'cm2'``.

    Returns
    -------
    factor : float
        The conversion factor.
    '''
    # To W=Joule/s => Joule = 1e7 erg
    factor = 1
    if m2_or_cm2 == 'cm2':
        factor = factor * 100 * 100
    factor = factor * 1e7  # erg/m2/s
    factor = factor / 1e4  # erg/cm2/s
    if sterad == False:
        beamsr = 1.133 * (resolution * const.a2r) ** 2
        factor = factor / beamsr  # erg/cm2/s/sr
    if nu_or_lambda == 'lambda':
        x = c/x
    # Umrechung zwischen ergs/s/cm2/sr = 2 k(CGS) nu^3/c(sm)^3 K km/s
    # => to make the units fit we have to multiply by 1*km in cm -> 1e5
    # i.e. const.km_in_cm
    # converts from K - > ergs
    conversionFactor = (2 * const.k_CGS * x ** 3 * const.km_in_cm /
                        (const.c_in_cm ** 3))
    factor = factor / conversionFactor
    if ToKKms == True:
        return 1/factor
    if ToKKms == False:
        return factor


def ergToKkms(x, toErg=False, nu_or_lambda='nu'):
    r"""
    Conversion between ergs/cm2/s/sr and K km/s.

    Parameters
    ----------
    x : float
        wavelenght/frequency [GHZ],
    toErg : True or False
        True converts the other direction, i.e. from K km/s to ergs/cm2/s/sr.
    nu_or_lambda : string
         Choose type of x: frequency = ``'nu'`` or wavelenght = ``'lambda'``.

    Returns
    -------
    factor : float
        The conversion factor.

    Notes
    -----
    Approved.
    """
    # To W=Joule/s => Joule = 1e7 erg
    factor = 1
    #print value
    if nu_or_lambda == 'lambda':
        x = const.c / x
    # Conversion between erg/s/cm2/sr = 2k(CGS) nu^3/c(cm)^3 K km/s
    # k(CGS) is Boltzsmanns constant in units of the CGS, nu the frequency of
    # the measusrement
    # c(cm) is the speed of light in cm.
    # => to make the units fit we have to multiply by 1*km in cm -> 1e5
    # i.e. const.km_in_cm
    # converts from K - > ergs
    conversionFactor = (2 * const.k_CGS * x ** 3 * const.km_in_cm /
                        (const.c_in_cm **3))
    factor = factor / conversionFactor
    if toErg == False:
        return factor
    if toErg == True:
        return 1 / factor


def Int2Lum(distance_in_pc, cm_or_m='cm'):
    r"""
    Conversion factor to calculate luminosity from intensities
    by integrating over the sky 4 pi Distance^2.

    Parameters
    ----------
    distance_in_pc : float
        Distance to the source in parsecs.
    cm_or_m : string
        Choose wether the out put is in cm^2 = ``'cm'`` or in
        m^2 = ``'m'``.

    Notes
    -----
    Approved.
    """
    if cm_or_m == 'm':
        return 4 * math.pi * (distance_in_pc * const.parsec_in_m_1) ** 2
    if cm_or_m == 'cm':
        return 4 * math.pi * ( distance_in_pc * const.parsec_in_cm) ** 2


def JyBToErgsB(input_flux, distance, wavelength, invert=False, map_use=False):
    r"""
    Conversion between Jy/beam and ergs/beam.

    Parameters
    ----------
    input_flux : float
        Flux to be converted in Jy/beam
    distance : float
        Distance to the source in parsec.
    wavelength : float
        Wavelength :math:`\lambda` in :math:`\mu m`.
    map_use :

    Returns
    -------
    The conversion factor (map_use = true) or the already converted flux
    (map_use = False).
    r"""
    # change from Jansky to erg s-1 cm-2 Hz-1
    conversion = JanskyToErgs_cm
    # integrate over sky ergs s-1 Hz-1
    conversion = conversion * Int2Lum(distance, cm_or_m='cm')
    # multiply by frequency
    conversion = conversion * const.c / (wavelength * 1e-6)
    if invert == False:
        if map_use == False:
            return input_flux * conversion
        if map_use == True:
            return conversion
    if invert == True:
        if map_use == False:
            return input_flux / conversion
        if map_use == True:
            return 1/conversion


def JyBToWM2Kpc2(input_flux, distance, major, minor, wavelength,
                 invert=False, map_use=False):
    r"""
    Conversion between Jy/beam and W m^-2 kpc^-2

    Parameters
    ----------
    input_flux :  float
        Flux to be converted.
    distance : float
        Distance to source in parsec.
    major : float
        Major Axis Beam (arcsec).
    minor : float
        Minor Axis Beam(arcsec).
    wavelength : float
        Wavelenght :math:`\lambda` in :math:`\mu m`
    invert : True or False
        Changes the direction of conversion.

    Returns
    -------
    float : the converted Flux.
    """
    # change to W/m2/Hz/beam
    conversion =  JanskyToWatt
    # calculate the beamsize in kpc2
    beamsize = 1.133 * (distance / 1e3) ** 2 * major * minor * const.a2r ** 2
    beamsInKpc2 = 1 / beamsize
    conversion = conversion * beamsInKpc2  # change to W/m2/Hz/kpc2
    conversion = conversion * const.c / (wavelength * 1e-6)  # to W/m2/kpc2
    if invert == False:
        return input_flux * conversion
    if invert == True:
        return input_flux / conversion


def JyBToWKpc2(input_flux, distance, major, minor,
               wavelength, invert=False, map_use=False):
    r"""
    Conversion from JyB to W kpc^-2.

    Parameters
    ----------
    input_flux :  float
        Flux to be converted.
    distance : float
        Distance to source in parsec.
    major : float
        Major Axis Beam (arcsec).
    minor : float
        Minor Axis Beam(arcsec).
    wavelength : float
        Wavelenght :math:`\lambda` in :math:`\mu m`.
    invert : True or False
        Changes the direction of conversion.

    Returns
    -------
    float : the converted Flux.
    """
    conversion =  JanskyToWatt  # change to W/m2/Hz/beam
    beamsize = (1.133 *
                (distance / 1e3) ** 2 *
                major *
                minor *
                const.a2r ** 2)  # calculate the beamsize in kpc2.
    beamsInKpc2 = 1 / beamsize
    conversion = conversion * beamsInKpc2  # change to W/m2/Hz/kpc2
    conversion = conversion * const.c / (wavelength * 1e-6)  # change to
                                                             # W/m2/kpc2
    conversion = conversion * Int2Lum(distance, cm_or_m='m')  # change to
                                                              # W/kpc2
    if not invert:
        return input_flux * conversion
    if invert:
        return input_flux / conversion

# The following functions are redundant maps/main.py contains the
# function flux_conversion, with the same functionality but with more 
# flexibility. This is however not usable outside of the maps environment.
# The following functions can be used generally and are correct!!

def kelvin_to_jansky(x, major, minor, nu_or_lambda='nu'):
    """
    Conversion from K.km/s (Tmb) and Jy/beam.

    Parameters
    ----------
    x : float
        wavelength/frequency [GHZ],
    major : float
        Major Axis Beam (arcsec),
    minor : float
        Minor Axis Beam(arcsec),
    nu_or_lambda : string
         Choose type of x: frequency = ``'nu'`` or wavelength = ``'lambda'``.

    Notes
    -----
    This function has been compared with the Time estimator from the
    [GILDAS] package ASTRO and yields the same conversion factors.

    References
    ----------
    .. [GILDAS] www.iram.fr/IRAMFR/GILDAS
    """
    if nu_or_lambda == 'lambda':
        def fcon(wavelengths, major, minor):
            return 1 / (1.359918e7 * wavelengths ** 2 / major / minor)
    if nu_or_lambda == 'nu':
        def fcon(frequency, major, minor):
            return 1 / (1.222233e6 * frequency ** (-2) / major / minor)
    return fcon(x, major, minor)


def jansky_to_kelvin(x, major, minor, nu_or_lambda='nu'):
    """
    Conversion from Jy/beam to K.km/s (Tmb).

    Parameters
    ----------
    x : float
        wavelength/frequency [GHZ],
    major : float
        Major Axis Beam (arcsec).
    minor : float
        Minor Axis Beam(arcsec).
    nu_or_lambda : string
         Choose type of x: frequency = ``'nu'`` or wavelength = ``'lambda'``.

    Notes
    -----
     Same as :func:`kelvin_to_jansky`
    """
    if nu_or_lambda == 'lambda':
        def fcon(wavelengths, major, minor):
            return 1.359918e7 * wavelengths **2 / major / minor
    if nu_or_lambda == 'nu':
        def fcon(frequency, Maj, Min):
            return 1.222233e6 * frequency ** (-2) / major / minor
    return fcon(x, major, minor)
