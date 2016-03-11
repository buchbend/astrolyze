from astropy import units as u
from astropy import constants as c

def Wm_to_Kkms(frequency=None, wavelenght=None, resolution=None):
    """Convert W/m^2 to K Km/s

    Conversion formula:

    .. math::

         1\,\frac{ergs}{s cm^2 sr} = \frac{2 k(CGS) nu^3}{c(sm)^3} K\,\frac{km}{s}

    """

    value = 1 * u.watt / u.meter**2
    if not frequency and not wavelenght:
        print "Please set frequency **or** wavelenght"
        raise SystemExit
    if frequency and wavelenght:
        print "Please set either frequency **or** wavelenght"
        raise SystemExit
    if frequency:
        frequency = frequency * u.GHz
    if wavelenght:
        wavelenght = wavelenght * u.meter
        frequency = wavelenght.to(u.GHz, equivalencies=u.spectral())
    if resolution == None:
        print "resolution needed"
        raise SystemExit

    # Convert W/m2 to erg/s/m2
    value = value.to(u.erg/u.meter**2/u.second)

    # Calculate the beamsize in sterad
    resolution = resolution * u.arcsec
    beamsr = 1.133 * (resolution.to(u.rad)) ** 2

    # 
    value = value / beamsr
    km_in_cm = (1 * u.km).to(u.cm)
    print km_in_cm
    c_in_cm = (c.c).to(u.cm/u.second)
    #=> to make the units fit we have to multiply by 1*km in cm -> 1e5
    #i.e. km_in_cm
    conversionFactor = (2 * c.k_B.cgs * frequency.to(u.Hz) ** 3 * km_in_cm /
                        (c_in_cm ** 3))
    value = value / conversionFactor
    return value.decompose() 

if __name__ == "__main__":
    conv = Wm_to_Kkms(wavelenght=122e-6, resolution=16)
    print 2.5 * conv
