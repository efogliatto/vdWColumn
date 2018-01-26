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








def vaporPhaseProfile(ci, Et, Eb, de, Tt, Tb):

    """
    Reduced density profile

    Arguments
    ci: initial condition. Concentration at Er_min
    Et: upper integration limit
    Eb: lower integration limit
    Tt: Temperature at Et
    Tb: Temperature at Eb
    de: energy step (estimated)
    """


    
    # Initial conditions

    nn = int(  np.floor( (Et - Eb) / de ) + 2  )

    Cr = ci * np.ones( nn )
    
    Er = Eb * np.ones( nn )

    dde = (Et - Eb) / (nn - 1)

    nablaTr = (Tt - Tb) / (Et - Eb)
    

    

    # Integrate differential equation

    for i in range(1,nn):


        Tr = Tb + nablaTr * (Er[i-1] - Eb)
        
        a = Tr / ((1.0 - Cr[i-1]/3.0)**2)

        b = 9.0 * Cr[i-1] / 4.0

        c = nablaTr * ( Cr[i-1] / (1 - Cr[i-1]/3.0) )

        
        f = ( Cr[i-1] + c )  /  ( a - b )
    
        Cr[i] = Cr[i-1] - dde * f

        Er[i] = Er[i-1] + dde
    




    # Integrate density profile

    mass = 0

    for i in range(1,nn):

        mass = mass + 0.5 * (Cr[i] + Cr[i-1])  * (Er[i] - Er[i-1])
    


    
    return Er, Cr, mass










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










    



# Reduced concentration profile

def rhoUniformTr( Tr = 0.99, c_bar = 1.0, Eb = 0., Et = 1., de = 1e-3 ):

    """
    Reduced concentration profile

    Arguments
    Eb: lower integration limit
    Et: upper integration limit
    de: energy step (estimated)
    """


    # Interphase densities

    cl, cg = interphaseDensities( Tr )


    # Initial mass excess

    mass = 1.


    # Initial positions

    Ett = Et

    Ebb = Eb

    Ei = 0.5 * Et  +  0.5 * Eb


    

    # Repeat until mass convergence

    while( abs(mass) > 1e-10 ):


        Er_g, C_g, mass_g = vaporPhaseProfile(cg, Et, Ei, de, Tr, Tr)

        Er_l, C_l, mass_l = liquidPhaseProfile(cl, Ei, Eb, de, Tr, Tr)


        mass = mass_g  +  mass_l  -  c_bar * (Et - Eb)


        # Too much liquid
    
        if ( mass > 0 ):

            Ett = Ei

            Ei = 0.5 * (Ei + Ebb)
            


        # Too much gas
        
        else:

            Ebb = Ei

            Ei = 0.5 * (Ett + Ei)            



    

    return Er_g, C_g, Er_l, C_l

    

    
