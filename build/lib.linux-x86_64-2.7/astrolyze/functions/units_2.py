import math

import astrolyze.functions.constants as const

# Constant conversion factors.
#==============> Approved !!! <==========================
WattToErgs    = 1e7  # 1W = 1e7 erg/s
ErgsToWatt    = 1e-7  # 1W = 1e-7 erg/s
JanskyToWatt  = 1e-26  # 1Jy = 1e-26 W/m2/Hz
WattToJansky  = 1e26  # 1W/m2/Hz  = 1e26 Jy
ErgsToJansky_cm  = 1e23  # 1 erg/s =  1e23 Jy * cm2 * Hz * s
JanskyToErgs_cm  = 1e-23  # 1 Jy = 1e-23 erg/s/cm2/Hz
ErgsToJansky_m  = 1e19  # 1 erg/s = 1e19 Jy * m2 * Hz * s
JanskyToErgs_m  = 1e-19  # 1 Jy = 1e-19 erg/s/m2/Hz
#==============> Approved !!! <==========================


# Natural Constants
c = 299792458.  # Speed of light [m]
c_in_cm = c * 1e2  # Speed of light [cm]
h = 6.62606896e-34  # Plancks constant [Js]
k = 1.3806503e-23  # Boltzman constant [m^2 kg s^-1 K^-1]
tBG = 2.7  # Cosmic Microwave Background Temperature in [K]
e = 2.7182818284  # Eulers number

# Natural Constants in cgs units.
k_CGS = 1.3806503e-16  # Boltzman constant [cm^2 g s^-1 K^-1]
h_CGS = 6.62606896e-27  # Plancks constant [Js]
c_CGS = 2.99792458e10  # Speed of light [cm]

# Angle Conversions
a2r = 4.848e-6  # arcsec to radian
a2r = 4.848e-6  # arcsec to radian
a2d = 1./60/60  # arcsec to degree
r2d = 180./math.pi  # radian to degree
r2a = 1./4.848e-6  # radian to arcsec
d2r = math.pi/180.  # degree to radian
d2a = 60*60  # degree to arcsec

def beam_in_sterad(resolution):
    beamsize = 1.133 * (a2r * resolution)**2
    print 'beam',beamsize
    return beamsize


def kelvin_to_jansky(x, major, minor, nu_or_lambda='nu'):
    """
    Conversion from K.km/s (Tmb) to Jy/beam.

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


def WmsrToKkms(frequency, beam):
    # W/m2/sr -> erg/s/m2/sr
    conversion = WattToErgs
    # erg/s/m2/sr -> erg/s/cm2/sr
    conversion = conversion / 100**2

    # Jy/beam -> K km/s
    conversion = conversion/kelvin_to_jansky(frequency/1e9,beam,beam)
    return conversion

print 0.6e-7*WmsrToKkms(1.4e12, 12)/100
print 1.2e-7*WmsrToKkms(1.4e12, 12)/100
print 2.1e-7*WmsrToKkms(1.4e12, 12)/100

