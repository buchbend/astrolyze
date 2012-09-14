# Copyright (C) 2012, Christof Buchbender
# BSD Licencse
import math

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

# Conversion of distances
parsec_in_m_1 = 3.08568025e16
parsec_in_m = 3.085e16  # parsec in m
parsec_in_cm = 3.08568025e18  # parsec in cm
km_in_cm = 1e5

# Masses
m_sun = 1.9891e30  # [kg]
m_proton = 1.672621637e-27  # [kg]

# Gauss constants
# GaussArea/(height*FWHM)
gauss_constant = 1.064467

# Luminosities 
LsunW = 3.846e26  # Watts
Lsunergs = 3.846e26*1e7  # erg/s
debye_to_EsuCm = 1.e-18  # Change from debye to esu/cm

# Angle Conversions
a2r = 4.848e-6  # arcsec to radian
a2d = 1./60/60  # arcsec to degree
r2d = 180./math.pi  # radian to degree
r2a = 1./4.848e-6  # radian to arcsec
d2r = math.pi/180.  # degree to radian
d2a = 60*60  # degree to arcsec

#################################################################
# Redundant but maybe still be used some in parts of the program.
pcInM = 3.085e16  # parsec in m
pcInCm = 3.08568025e18  # parsec in cm
arcsecInRad = 4.848e-6
arcsecInGrad = 1./60./60
squareArcsecInSterad = 4.254517e10
gaussConst = 1.064467
mSun = 1.9891e30  # [kg]
mProton = 1.672621637e-27  # [kg]
parsecInMeter = 3.08568025e16  # parsec in m
parsecInCentiMeter = 3.08568025e18  # parsec in cm
