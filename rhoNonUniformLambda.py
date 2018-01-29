from math import log

import numpy as np

from scipy.optimize import fsolve

from scipy.interpolate import interp1d

from numpy.linalg import norm, solve

from .interphaseDensities import interphaseDensities

from .vaporPhaseProfileWithT import vaporPhaseProfileWithT

from .liquidPhaseProfileWithT import liquidPhaseProfileWithT






def updateTemperature(Er, Tr, Cl, Cg, idx, kappa = 1.0):

    """
    Update Temperature distribution

    Arguments
    Er: reduced height grid
    Tr: initial Temperature distribution
    Cl: liquid density
    Cg: vapor density
    """


    def lambdaRho(rho,kappa):

        return kappa * rho

        

    dde = Er[1] - Er[0]

    
    
    # Liquid phase matrix

    r = lambdaRho( Cg[idx], kappa ) * (Tr[idx+1] - Tr[idx]) / dde
    
    Ml = np.zeros( (idx+1, idx+1) )

    Bl = r * np.ones(idx+1)
    



    # vapor phase matrix

    r = lambdaRho( Cl[idx], kappa ) * (Tr[idx] - Tr[idx-1]) / dde    

    Mg = np.zeros( (len(Er)-idx, len(Er)-idx) )

    Bg = r * np.ones( len(Er)-idx )

    
    beta = 1 / (lambdaRho(Cl[idx],kappa) - lambdaRho(Cg[idx],kappa) )



    Told = np.ones( len(Tr) )

    err = 1.


    # while( np.linalg.norm(Tr - Told) > 1e-3 ):
    for k in range(10):

        
        # Liquid
        
        for i in range(idx+1):

        
            if i == 0:

                Ml[i,i] = -lambdaRho( Cl[i], kappa ) / dde

                Ml[i,1] = lambdaRho( Cl[i], kappa ) / dde



            elif i == idx:

                Ml[i,i] = 1.

                Ml[i,i-1] = -beta * lambdaRho(Cl[i],kappa)

                Bl[i] = beta * lambdaRho(Cg[i],kappa) * Tr[idx+1]



            else:

                Ml[i,i+1] = lambdaRho( Cl[i], kappa ) / (2*dde)

                Ml[i,i-1] = -lambdaRho( Cl[i], kappa ) / (2*dde)




        Tl = solve(Ml, Bl)

        for i in range(idx+1):

            Tr[i] = Tl[i]


        


        # Vapor
            

        for i in range( len(Er)-idx ):

        
            if i == 0:

                Mg[i,i] = 1.

                Mg[i,i+1] = -beta * lambdaRho(Cg[i+idx],kappa)

                Bg[i] = beta * lambdaRho(Cl[i+idx],kappa) * Tr[idx-1]
            

            
            elif i == len(Er) - idx -1:

                Mg[i,i] = lambdaRho( Cg[i+idx], kappa ) / dde

                Mg[i,i-1] = -lambdaRho( Cg[i+idx], kappa ) / dde



            else:

                Mg[i,i+1] = lambdaRho( Cg[i+idx], kappa ) / (2*dde)

                Mg[i,i-1] = -lambdaRho( Cg[i+idx], kappa ) / (2*dde)
            


            
        Tg = solve(Mg, Bg)

        for i in range( len(Er)-idx ):

            Tr[i+idx] = Tg[i]

            

            
        err = np.linalg.norm(Tr - Told)

        for i in range(len(Tr)):

            Told[i] = Tr[i]

            
        print(err)

    

    pass
    



    



 

# Reduced concentration profile. Fixed Gradient

def rhoNonUniformLambda( Tt = 0.99, Tb = 0.99, c_bar = 1.0, Eb = 0., Et = 1e-03, npoints = 10000, kappa = 1.0 ):

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

    # for k in range(1):
    
    while( abs(mass) > ((Et-Eb)/npoints) ):


        # Interphase densities

        cl, cg = interphaseDensities( Tr[Ei] )
        
        
        # Liquid an vapor profiles

        C_g, mass_g = vaporPhaseProfileWithT(cg, Er, Tr, Ei)

        C_l, mass_l = liquidPhaseProfileWithT(cl, Er, Tr, Ei)

        
        
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



            
    # Update temperature distribution

    updateTemperature(Er, Tr, C_l, C_g, Ei)        

    

    return Er, C_g, C_l, Ei




