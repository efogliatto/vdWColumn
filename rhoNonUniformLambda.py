from math import log

import numpy as np

from scipy.optimize import fsolve

from scipy.interpolate import interp1d

from numpy.linalg import norm, solve

from .interphaseDensities import interphaseDensities

from .vaporPhaseProfileWithT import vaporPhaseProfileWithT

from .liquidPhaseProfileWithT import liquidPhaseProfileWithT
    
from .updateTemperature import updateTemperature


 

# Reduced concentration profile. Fixed Gradient

def rhoNonUniformLambda( Tt = 0.99,
                         Tb = 0.99,
                         c_bar = 1.0,
                         Eb = 0.,
                         Et = 1e-03,
                         npoints = 10000,
                         kappa = 1.0,
                         updateT = False,
                         thcond = 'uniform',
                         ttol = 1e-3):

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

    Tnew = np.linspace(Tb, Tt, npoints)

    Terr = 1.


    
    # Initial positions (as index)

    Ett = npoints

    Ebb = 0

    Ei = int( npoints/2 )




    

    # Repeat until mass convergence

    # for k in range(1):
    
    while( abs(mass) > ((Et-Eb)/npoints) ):


        Terr = 1.0
        

        # First, iterate over c and T distributions

        while( Terr > 1e-3 ) :

            
            # Interphase densities

            cl, cg = interphaseDensities( Tr[Ei] )

            
        
            # Liquid an vapor profiles

            C_g, mass_g = vaporPhaseProfileWithT(cg, Er, Tr, Ei)

            C_l, mass_l = liquidPhaseProfileWithT(cl, Er, Tr, Ei)

            

            # Update temperature distribution

            if updateT == True:
    
                Tnew = updateTemperature(Er, C_l, C_g, Ei, Tt, Tb, kappa, thcond)            

                Terr = np.linalg.norm(Tnew - Tr)
            

                for i in range(npoints):

                    Tr[i] = Tnew[i]

            else:

                Terr = 0.





        
        
        # Mass excess
        
        mass = mass_g  +  mass_l  -  c_bar * (Et - Eb)

        

        # Too much liquid
    
        if ( mass > 0 ):

            Ett = Ei

            Ei = int( 0.5 * (Ei + Ebb) )
            


        # Too much gas
        
        else:

            Ebb = Ei

            Ei = int( 0.5 * (Ett + Ei) )




    return Er, C_g, C_l, Ei, Tr
