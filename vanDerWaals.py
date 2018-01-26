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




# Reduced concentration profile

def cprofile( ci, Tr = 0.99, Er_min = 0., Er_max = 1., de = 1e-3 ):

    """
    Reduced concentration profile

    Arguments
    ci: initial condition. Concentration at interphase
    Er_min: interphase position
    Er_max: integration limit
    de: energy step (estimated)
    """


    # Initial conditions

    nn = int(  np.floor( (Er_max - Er_min) / de ) + 2  )

    C = ci * np.ones( nn )

    Er = np.zeros( nn )



   
    # M construction

    M = np.zeros( (nn, nn) )
    
    for i in range(nn):

        if i == 0:

            M[i][i]    =  1


        elif i == (nn-1):

            M[i][i]    =  3./2.

            M[i][i-1] =  -4./2.

            M[i][i-2] =   1./2.


        else:

            M[i][i+1]  =   1./2.

            M[i][i-1]  =  -1./2.
    






    # Equation function

    def fc( A, Tr, c ):
    
        CC = - A / ( (Tr / (1-A/3.)**2) - (9./4.) * A  )

        CC[0] = c

        return CC



    # Solve equation

    C_new = np.zeros( nn )
    
    b = de * fc( C, Tr, ci )
    
    
    
    while( np.fabs( norm( C_new - C ) )   > 1e-1 ):

        C = C_new
        
        C_new = solve(M, b)

        b = fc( C_new, Tr, ci )

        



        
    
            


    pass

    

    
