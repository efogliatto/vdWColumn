from math import log

import numpy as np

from scipy.optimize import fsolve

from scipy.interpolate import interp1d

from numpy.linalg import norm, solve

from .interphaseDensities import interphaseDensities

from .vaporPhaseProfileWithT import vaporPhaseProfileWithT

from .liquidPhaseProfileWithT import liquidPhaseProfileWithT






def updateTemperature(Er, Tr, Cl, Cg, idx, kappa = 1.0, thcond = 'uniform', ttol = 1e-3):

    """
    Update Temperature distribution

    Arguments
    Er: reduced height grid
    Tr: initial Temperature distribution
    Cl: liquid density
    Cg: vapor density
    """


    # Thermal conductivity
    
    def lambdaRho( rho, kappa, thcond = 'uniform' ):

        l = kappa

        
        if thcond == 'uniform':

            l = kappa

        elif thcond == 'linear':

            l = kappa * rho

            
        return l




    # Reduced energy step

    dde = Er[1] - Er[0]

    
    
    # Liquid phase matrix
    
    r = lambdaRho( Cg[idx], kappa, thcond ) * (Tr[idx+1] - Tr[idx]) / dde
    
    Ml = np.zeros( (idx+1, idx+1) )

    Bl = np.ones(idx+1)
    



    # vapor phase matrix

    r = lambdaRho( Cl[idx], kappa, thcond ) * (Tr[idx] - Tr[idx-1]) / dde    

    Mg = np.zeros( (len(Er)-idx, len(Er)-idx) )

    Bg = np.ones( len(Er)-idx )

       



    Told = np.ones( len(Tr) )

    err = 1.


    while( err > ttol ):
    # for k in range(10):


        beta = 1 / (lambdaRho(Cl[idx], kappa, thcond) + lambdaRho(Cg[idx], kappa, thcond) )
        

        # Liquid
        
        
        # Flux
        
        r = lambdaRho( Cg[idx], kappa, thcond ) * (Tr[idx+1] - Tr[idx]) / dde
        
        Bl.fill( r )
        
        
        # Matrix
        
        for i in range(idx+1):

        
            if i == 0:

                Ml[i,i] = -lambdaRho( Cl[i], kappa, thcond ) / dde

                Ml[i,1] = lambdaRho( Cl[i], kappa, thcond ) / dde



            elif i == idx:

                Ml[i,i] = 1.

                Ml[i,i-1] = -beta * lambdaRho(Cl[i], kappa, thcond)

                Bl[i] = beta * lambdaRho(Cg[i], kappa, thcond) * Tr[idx+1]



            else:

                Ml[i,i+1] = lambdaRho( Cl[i], kappa, thcond ) / (2*dde)

                Ml[i,i-1] = -lambdaRho( Cl[i], kappa, thcond ) / (2*dde)




        Tl = solve(Ml, Bl)

        for i in range(idx+1):

            Tr[i] = Tl[i]


        


        # Vapor


        # Flux

        r = lambdaRho( Cl[idx], kappa, thcond ) * (Tr[idx] - Tr[idx-1]) / dde    

        Bg.fill(r)


        
        # Matrix

        for i in range( len(Er)-idx ):

        
            if i == 0:

                Mg[i,i] = 1.

                Mg[i,i+1] = -beta * lambdaRho(Cg[i+idx], kappa, thcond)

                Bg[i] = beta * lambdaRho(Cl[i+idx], kappa, thcond) * Tr[idx-1]
            

            
            elif i == len(Er) - idx -1:

                Mg[i,i] = lambdaRho( Cg[i+idx], kappa, thcond ) / dde

                Mg[i,i-1] = -lambdaRho( Cg[i+idx], kappa, thcond ) / dde



            else:

                Mg[i,i+1] = lambdaRho( Cg[i+idx], kappa, thcond ) / (2*dde)

                Mg[i,i-1] = -lambdaRho( Cg[i+idx], kappa, thcond ) / (2*dde)
            


            
        Tg = solve(Mg, Bg)

        for i in range( len(Er)-idx ):

            Tr[i+idx] = Tg[i]

            


        # Compute error
            
        # err = np.linalg.norm(Tr - Told)

        err = np.linalg.norm(Tr - Told) / np.linalg.norm(Tr)

        

        for i in range(len(Tr)):

            Told[i] = Tr[i]
            

        # print(err)
    

    return Told
    



    



 

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

        if updateT == True:
    
            Tr = updateTemperature(Er, Tr, C_l, C_g, Ei, kappa, thcond, ttol)


            

    return Er, C_g, C_l, Ei, Tr




