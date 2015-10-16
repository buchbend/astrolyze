import math

def calc_absolute_magnitude(apparent_magnitude, distance):
    """ calculate absolute mag. from apparent mag. and distance"""
    m = apparent_magnitude
    d = distance
    return m - 5 * math.log(d,10) + 5

def lj_to_pc(lj):
    """Convert lightyears to parsec"""
    return lj/3.26
