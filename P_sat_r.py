from math import log

import numpy as np

from scipy.optimize import fsolve

from scipy.interpolate import interp1d

from numpy.linalg import norm, solve



# Reduced saturation pressure

def P_sat_r ( cg, cl ):

    """Saturation pressure
       
       Arguments
       cg: reduced concentration of vapor phase
       cl: reduced concentration of liquid phase"""
    
    return cg * cl * (3.0 - (cg + cl))
