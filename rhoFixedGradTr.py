from math import log

import numpy as np

from scipy.optimize import fsolve

from scipy.interpolate import interp1d

from numpy.linalg import norm, solve

from interphaseDensities import interphaseDensities

from vaporPhaseProfile import vaporPhaseProfile

from liquidPhaseProfile import liquidPhaseProfile




# Reduced concentration profile. Fixed Gradient

def rhoFixedGradTr( Tt = 0.99, Tb = 0.99, c_bar = 1.0, Eb = 0., Et = 1., de = 1e-3 ):

    """
    Reduced concentration profile

    Arguments
    Eb: lower integration limit
    Et: upper integration limit
    de: energy step (estimated)
    """


    # Interphase densities

    cl, cg = interphaseDensities( 0.5*Tt  +  0.5*Tb )


    # Initial mass excess

    mass = 1.


    # Initial positions

    Ett = Et

    Ebb = Eb

    Ei = 0.5 * Et  +  0.5 * Eb



    nablaTr = (Tt - Tb) / (Et - Eb)
    

    # Repeat until mass convergence

    while( abs(mass) > 1e-10 ):


        Er_g, C_g, mass_g = vaporPhaseProfile(cg, Et, Ei, de, Tt, Tb + nablaTr*(Ei-Eb))

        Er_l, C_l, mass_l = liquidPhaseProfile(cl, Ei, Eb, de, Tb + nablaTr*(Ei-Eb), Tb)


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
