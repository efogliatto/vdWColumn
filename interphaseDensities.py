from math import log

import numpy as np

from scipy.optimize import fsolve

from scipy.interpolate import interp1d

from numpy.linalg import norm, solve

from .cequation import cequation


# Reduced concentrations for specific reduced temperature

def interphaseDensities( Tr ):

    """
    Reduced concentration, both phases.

    Arguments
    Tr: reduced temperature
    """


    # Initial conditions. Helps convergence

    itr = np.array([0.99, 0.95, 0.93, 0.90, 0.70, 0.60, 0.50, 0.45, 0.40])

    icl = np.array([1.1, 1.46, 1.54, 1.65, 1.9, 2.2, 2.4, 2.5, 2.7])

    icg = np.array([0.8, 0.58, 0.51, 0.43, 0.1, 0.05, 0.025, 0.0125, 1.5625e-03])

    

    # Solve cequation to find concentrations

    fl = interp1d(itr, icl)

    fg = interp1d(itr, icg)
    
    cl, cg = fsolve( cequation, [fl(Tr), fg(Tr)], (Tr, ), xtol=1e-10 )
    

    
    return cl, cg
