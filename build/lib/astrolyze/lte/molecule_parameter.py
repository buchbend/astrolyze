# Copyright (C) 2012, Christof Buchbender
# BSD Licencse
r"""
This script generates a dictionary storing the information of molecular
transitions.
"""
from numpy import interp

import astrolyze.functions.constants as const

class Molecule:
    r"""
    A class holding attributes that define the characteristics of an
    individual tansition of a molecules.

    The :py:func:`astrolyze.functions.astro_functions.calc_N` routine
    depends on this class.

    Parameters
    ----------

    nu : float
    Q : list
        The partition function of the molecule evaluated at discreet
        temperatures. Given in T.
    T : list
        The temperatures at which the partition functions was evaluated. The
        value of Q used finally is interpolated to the excitation_temperature
        give.
    Eu : float
        Energy of the upper state of the transition.
    Aul : float
        The Einstein coefficient. For the transition of upper-to-lower state.
    gu : float
        ADD DESCRIPTION.
    mu : float
        ADD DESCRIPTION.
    name : string
        The name of the molecule. 
    """
    def __init__(self, nu, Q, T, Eu, Aul, gu, mu, name=''):
        self.nu = nu 
        self.Q = Q
        self.T = T
        self.Eu = Eu
        self.Aul = Aul
        self.gu = gu
        self.name = name
        self.mu = mu

molecule_dictionary = {}

# HCO+ (1-0)
nu  =  89.1885247e9   # Hz
# Z the partition function has to be interpolated to the desired value of T 
Q = [2.4426, 2.1638, 2.0276, 1.8481, 1.5488, 1.2519, 0.9592, 0.6748, 0.4315,
     0.2225]
T   = [500.0, 300.0, 225.0, 150.0, 75.0, 37.50, 18.75, 9.375, 5.0, 2.725]
Eu  = 2.9750
Aul =  10 ** (-2.2608)
gu  = 3.
mu  = 3.9*const.debye_to_EsuCm
HCO = Molecule(nu, Q, T, Eu, Aul, gu, mu, 'HCO+(1-0)')
molecule_dictionary['HCO'] = HCO 

# HOC+ (1-0)
nu  =  89.4874140e9  # Hz
# Z the partition function has to be interpolated to the desired value of T 
Q   = [2.3558, 2.1739, 1.9228, 1.5554, 1.2505, 0.9578, 0.6735, 0.4303, 0.2214]
T   = [300.0, 225.0, 150.0, 75.0, 37.50, 18.75, 9.375, 5.0, 2.725]
Eu  = 2.9850
Aul = 10 ** (-2.7467)
gu  = 3.
mu  = 2.77*const.debye_to_EsuCm
HOC = Molecule(nu, Q, T, Eu, Aul, gu, mu, 'HOC+(1-0)')
molecule_dictionary['HOC']=HOC 

# HCN (1-0)
nu  = 88.6316022e9  # Hz
Q = [2.9688, 2.6566, 2.5122, 2.3286, 2.0286, 1.7317, 1.4389, 1.1545, 0.9109,
     0.7016]
T   = [500.0, 300.0, 225.0, 150.0, 75.0, 37.50, 18.75, 9.375, 5.0, 2.725]
Eu  = 2.9564
Aul = 10 ** (-2.5140)
gu  = 9.
mu  = 2.985*const.debye_to_EsuCm
HCN = Molecule(nu, Q, T, Eu, Aul, gu, mu, 'HCN(1-0)')
molecule_dictionary['HCN'] = HCN 

# HNC (1-0)
nu  = 90.6635680e9  # Hz
Q = [2.5455, 2.2255, 2.0585, 1.8507, 1.5419, 1.2449, 0.9523, 0.6683, 0.4254,
     0.2174]
T   = [500.0, 300.0, 225.0, 150.0, 75.0, 37.50, 18.75, 9.375, 5.0, 2.725]
Eu  = 3.0242
Aul = 10 ** (-2.5218)
gu  = 3.
mu  = 3.05*const.debye_to_EsuCm
HNC = Molecule(nu, Q, T, Eu, Aul, gu, mu, 'HNC(1-0)')
molecule_dictionary['HNC'] = HNC 

# 12CO (1-0)
nu  =  115.271e9  # Hz
Q = [2.2584, 2.0369, 1.9123, 1.7370, 1.4386, 1.1429, 0.8526, 0.5733, 0.3389,
     0.1478]
T   = [500.0, 300.0, 225.0, 150.0, 75.0, 37.50, 18.75, 9.375, 5.0, 2.725]
Eu  = 3.8450
Aul = 10 ** (-5.0105)
gu  = 3.
mu  = 0.11011*const.debye_to_EsuCm
CO1210 = Molecule(nu, Q, T, Eu, Aul, gu, mu, '12(1-0)')
molecule_dictionary['CO1210'] = CO1210 

# 13CO (1-0)
nu = 110.2013543e9   # Hz
Q = [2.5789, 2.3574, 2.2328, 2.0575, 1.7589, 1.4630, 1.1722, 0.8919, 0.6558,
     0.4611]
T = [500.0, 300.0, 225.0, 150.0, 75.0, 37.50, 18.75, 9.375, 5.0, 2.725]
Eu = 3.6759
Aul = 10 ** (-5.0662)
gu = 6
mu = 0.11046*const.debye_to_EsuCm
CO1310 = Molecule(nu, Q, T, Eu, Aul, gu, mu, '13CO(1-0)')
molecule_dictionary['CO1310'] = CO1310
molecule_dictionary['13CO'] = CO1310

# 12CO (2-1)
nu = 230.537990e9   # Hz
Q = [2.2584, 2.0369, 1.9123, 1.7370, 1.4386, 1.1429, 0.8526, 0.5733, 0.3389,
     0.1478]
T = [500.0, 300.0, 225.0, 150.0, 75.0, 37.50, 18.75, 9.375, 5.0, 2.725]
Eu = 11.5350
Aul = 10 ** (-4.1197)
gu = 5.
mu = 0.11011 * const.debye_to_EsuCm 
CO1221 = Molecule(nu, Q, T, Eu, Aul, gu, mu, '12CO(2-1)')
molecule_dictionary['CO1221'] = CO1221 
molecule_dictionary['12CO'] = CO1221 

# CN 
nu = 112.1016560e9   # Hz
Q = [3.0450, 2.8222, 2.6976, 2.5223, 2.2238, 1.9280, 1.6376, 1.3579, 1.1228,
     0.9303]
T = [500.0, 300.0, 225.0, 150.0, 75.0, 37.50, 18.75, 9.375, 5.0, 2.725]
Eu = 2042.4222
Aul = 10 ** (-8.0612)
gu = 2.
mu = 1.45 * const.debye_to_EsuCm 
CN = Molecule(nu, Q, T, Eu, Aul, gu, mu, '12CO(2-1)')
molecule_dictionary['CN'] = CN

# CCH 
# CCH has six hyperfine transitions, the values given here correspond to the
# strongest transition.
nu = 87.3168980e9  # Hz
Q = [3.9556, 3.3192, 2.9114, 2.7164, 2.4836, 2.1605, 1.8628, 1.5699, 1.2852,
     1.0411, 0.8308]
T = [1000.0, 500.0, 300.0, 225.0, 150.0, 75.0, 37.50, 18.75, 9.375, 5.0, 2.725]
Eu = 2.914
Aul = 10 ** (-5.2060)
gu = 5.
mu = 0.77 * const.debye_to_EsuCm 
CCH = Molecule(nu, Q, T, Eu, Aul, gu, mu, 'CCH(1-0)')
molecule_dictionary['CCH'] = CCH  






