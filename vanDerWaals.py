from math import log

import numpy as np

from scipy.optimize import fsolve

from scipy.interpolate import interp1d



# Reduced saturation pressure

def P_sat_r ( cg, cl ):

    """Saturation pressure
       
       Arguments
       cg: reduced concentration of vapor phase
       cl: reduced concentration of liquid phase"""
    
    return cg * cl * (3.0 - (cg + cl))




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




# Reduced concentrations for specific reduced temperature

def c_r( Tr ):

    """
    Reduced concentration, both phases.

    Arguments
    Tr: reduced temperature
    """


    # Initial conditions. Helps convergence

    itr = np.array([0.99, 0.95, 0.90, 0.70, 0.60, 0.50, 0.45, 0.40])

    icl = np.array([1.1, 1.3, 1.5, 1.9, 2.2, 2.4, 2.5, 2.7])

    icg = np.array([0.8, 0.5, 0.3, 0.1, 0.05, 0.025, 0.0125, 1.5625e-03])

    

    # Solve cequation to find concentrations

    fl = interp1d(itr, icl)

    fg = interp1d(itr, icg)
    
    cl, cg = fsolve( cequation, [fl(Tr), fg(Tr)], (Tr, ), xtol=1e-10 )
    

    
    return cl, cg
