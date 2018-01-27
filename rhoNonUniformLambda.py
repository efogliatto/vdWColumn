from math import log

import numpy as np

from scipy.optimize import fsolve

from scipy.interpolate import interp1d

from numpy.linalg import norm, solve

from interphaseDensities import interphaseDensities

from vaporPhaseProfileWithT import vaporPhaseProfileWithT

from liquidPhaseProfileWithT import liquidPhaseProfileWithT


 

# Reduced concentration profile. Fixed Gradient

def rhoNonUniformLambda( Tt = 0.99, Tb = 0.99, c_bar = 1.0, Eb = 0., Et = 1e-03, npoints = 10 ):

    """
    Reduced concentration profile

    Arguments
    Eb: lower integration limit
    Et: upper integration limit
    de: energy step (estimated)
    """

    

    # Initial mass excess

    mass = 1.


    # Initial distributions

    Er = np.linspace(Eb, Et, npoints)

    Tr = np.linspace(Tb, Tt, npoints)


    
    # Initial positions (as index)

    Ett = npoints

    Ebb = 0

    Ei = int( npoints/2 )




    

    # Repeat until mass convergence

    for k in range(1):
    
    # while( abs(mass) > 1e-10 ):


        # Interphase densities

        cl, cg = interphaseDensities( Tr[Ei] )
        print(cl)
        print(cg)
        
        
        # Liquid an vapor profiles

        C_g, mass_g = vaporPhaseProfileWithT(cg, Er, Tr, Ei)

        C_l, mass_l = liquidPhaseProfileWithT(cl, Er, Tr, Ei)

        
        
        # Mass excess
        
        mass = mass_g  +  mass_l  -  c_bar * (Et - Eb)

        # print(mass)

        # Too much liquid
    
        if ( mass > 0 ):

            Ett = Ei

            Ei = int( 0.5 * (Ei + Ebb) )
            


        # Too much gas
        
        else:

            Ebb = Ei

            Ei = int( 0.5 * (Ett + Ei) )




        # # New interphase temperature

        # hl = Ei - Eb

        # hg = Et - Ei
        
        # Ti = ( Tb / hl  +  gamma * Tt / hg)  *  ( hg * hl / ( gamma * hl + hg ) )
            

    

    return Er, C_g, C_l, Ei
