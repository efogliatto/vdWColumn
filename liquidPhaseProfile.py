from math import log

import numpy as np

from scipy.optimize import fsolve

from scipy.interpolate import interp1d

from numpy.linalg import norm, solve




def liquidPhaseProfile(ci, Et, Eb, de, Tt, Tb):

    """
    Reduced density profile

    Arguments
    ci: initial condition. Concentration at Et
    Et: upper integration limit
    Eb: lower integration limit
    Tt: Temperature at Et
    Tb: Temperature at Eb
    de: energy step (estimated)
    """


    
    # Initial conditions

    nn = int(  np.floor( (Et - Eb) / de ) + 2  )

    Cr = ci * np.ones( nn )
    
    Er = Et * np.ones( nn )

    dde = (Et - Eb) / (nn - 1)

    nablaTr = (Tt - Tb) / (Et - Eb)
    
    
    

    # Integrate differential equation

    for i in range(nn-2,-1,-1):


        Tr = Tb + nablaTr * (Er[i+1] - Eb)
        
        a = Tr / ((1.0 - Cr[i+1]/3.0)**2)

        b = 9.0 * Cr[i+1] / 4.0

        c = nablaTr * ( Cr[i+1] / (1 - Cr[i+1]/3.0) )

        
        f = ( Cr[i+1] + c )  /  ( a - b )
    
        Cr[i] = Cr[i+1] + dde * f

        Er[i] = Er[i+1] - dde

        if(Er[i] < 0):

            Er[i] = 0.
    




    # Integrate density profile

    mass = 0

    for i in range(1,nn):

        mass = mass + 0.5 * (Cr[i] + Cr[i-1]) * (Er[i] - Er[i-1])
    


    
    return Er, Cr, mass
