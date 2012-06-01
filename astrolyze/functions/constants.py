import math
# natural Constants

c = 299792458.	#Speed of light [m]
h = 6.62606896e-34 # Plancks constant [Js]
k = 1.3806503e-23 # Boltzman constant [m^2 kg s^-1 K^-1]
tBG = 2.7 # Cosmic Microwave Background Temperature in [K]

# in CGS
k_CGS = 1.3806503e-16 # Boltzman constant [cm^2 g s^-1 K^-1]
h_CGS = 6.62606896e-27 # Plancks constant [Js]
c_CGS=2.99792458e10  #Speed of light [cm]

e=2.7182818284 # Eulers number 

# Distances

parsecInMeter = 3.08568025e16 # parsec in m
parsecInCentiMeter = 3.08568025e18# parsec in cm

# redundant but maybe used in program parts.
pcInM = 3.085e16# parsec in m
pcInCm = 3.08568025e18 # parsec in cm

arcsecInRad = 4.848e-6 
arcsecInGrad = 1./60./60

squareArcsecInSterad = 4.254517e10

#Masses

mSun = 1.9891e30 # [kg]
mProton = 1.672621637e-27 #[kg]


#Gauss constants 
# GaussArea/(height*FWHM)
gaussConst = 1.064467

# Luminosities

LsunW = 3.846e26 # [W]
Lsunergs = 3.846e26*1e7 # erg/s 

debye_to_EsuCm = 1.e-18  # change from debye to esu/cm

# angle Conversions

# a: arseconds
# g: grad
# d: degrees
# r: radian

a2r = 4.848e-6
a2d = 1./60/60

r2d = 180./math.pi
r2a = 1./4.848e-6

d2r = math.pi/180.
d2a = 60*60








