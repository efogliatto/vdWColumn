from math import log

import numpy as np

from scipy.optimize import fsolve

from scipy.interpolate import interp1d

from numpy.linalg import norm, solve




def vaporPhaseProfileWithT(ci, Er, Tr, idx):

    """
    Reduced density profile on fixed grid

    Arguments
    ci: initial condition. Concentration at Er[idx]
    Er: reduced height grid
    Tr: Temperature distribution (on Er grid)
    """


    
    # Initial conditions

    Cr = ci * np.ones( len(Er) )

    dde = Er[1] - Er[0]    

    

    # Integrate differential equation only vaporPhase

    for i in range(idx + 1, len(Er)):
        
        
        a = Tr[i-1] / ((1.0 - Cr[i-1]/3.0)**2)

        b = 9.0 * Cr[i-1] / 4.0

        nablaTr = (Tr[i] - Tr[i-2]) / (2*dde)
        
        c = nablaTr * ( Cr[i-1] / (1 - Cr[i-1]/3.0) )

        
        f = ( Cr[i-1] + c )  /  ( a - b )
    
        Cr[i] = Cr[i-1] - dde * f
    




    # Integrate density profile

    mass = 0

    for i in range(idx + 1, len(Er)):

        mass = mass + 0.5 * (Cr[i] + Cr[i-1])  * (Er[i] - Er[i-1])
    


    
    return Cr, mass
