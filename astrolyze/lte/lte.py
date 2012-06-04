r"""
Functions to calculate LTE column densities.

Add Documentation.
"""
from numpy import interp
import astrolyze.functions.constants as const

def calc_jnu(nu, T):
    r"""
    Calculates :math:`J_{\nu}` needed for lte_column_density.
    !!! CGS Units!!! Extend Documentations!!!

    Parameters
    ----------

    nu: float
        Frequency 
    T: float
        Temperature
    """
    return (const.h_CGS*nu / const.k_CGS / (exp(const.h_CGS*nu / const.k_CGS /
            T) - 1))

def lte_column_density(nu, Tmb, excitation_temperature, J, Z, mu):
    """
    This function calculates the Column densities of linear molecules 

    Units are all to be given in cgs
    Z is the array of partition function values for the corresponding
    temperatures in T these are the log values of Z 

    Notes
    -----

    Extend documentation!!!!
    """
    print 'excitation_temperature', excitation_temperature
    print 'Tmb', Tmb
    hNuKT = (const.h_CGS * nu) / (const.k_CGS * excitation_temperature)
    colDens = 3 * const.h_CGS / (8 * math.pi**3 * mu**2)
    colDens *= Z / J
    colDens *= exp(hNuKT)
    colDens *= 1 / (1 - exp( - 1 * hNuKT))
    colDens *= 1 / (calc_jnu(nu,excitation_temperature) - Jnu(nu, const.tBG))
    colDens *= Tmb * 1e5
    return colDens

def calc_N(molecule, excitation_temperature, J, W):
    r"""
    Calculates the column density for a molecule. 
    !!! LOOK into the remaining Code and merge!!!
    """
    T=[]
    # reverse the Arrays T and Q of the molecules for the interpolating function
    for i in range(len(molecule.T)):
        T+=[molecule.T[len(molecule.T)-i-1]]
    Q=[]
    for i in range(len(molecule.Q)):
        Q+=[molecule.Q[len(molecule.Q)-i-1]]
    # interpolate the partition function for excitation_temperature from the values provided 
    # by CDMS   
    Z = interp(excitation_temperature, T, Q)
    print Z
    Z = 10**Z            # change from log Z to Z
    return lte_column_density(molecule.nu, W, excitation_temperature, J, Z, molecule.mu)

def calc_excitation_temperature(Tb, nu):
    """
    Calculation of the excitation temperature of an optically thick 12CO line 
    under the assumption of LTE.

    Parameters
    ----------

    Tb:
    """
    excitation_temperature = (const.h_CGS * nu)/const.k_CGS
    excitation_temperature *= ((log( 1 + (((const.k_CGS * Tb / const.h_CGS / nu) + (1 /
           (exp((const.h_CGS * nu)/(const.k_CGS*const.tBG)) - 1))) ** (-1))))
           ** (-1))
    return excitation_temperature