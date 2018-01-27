from math import log

import numpy as np

from scipy.optimize import fsolve

from scipy.interpolate import interp1d

from numpy.linalg import norm, solve

from P_sat_r import P_sat_r



# Reduced concentration function. Used for non-linear equation

def cequation( c, Tr ):

    """Reduced concentration. Evaluation of system of equations
       
       Arguments
       cg: reduced concentration of vapor phase
       cl: reduced concentration of liquid phase
       Tr: reduced temperature"""

    # Unpack values
    cl, cg = c

    # Maxwell rule
    y0 = 8.0 * Tr * log( (3.0/cg - 1.0) / (3.0/cl - 1.0) )   +   9.0 * (cg - cl)  -  3.0 * P_sat_r(cg,cl) * (1.0/cg - 1.0/cl)

    # vdW equation
    y1 = 8.0 * Tr - (cg+cl) * (3.0-cl) * (3.0-cg)

  
    return y0, y1
