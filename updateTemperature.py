import numpy as np

from numpy.linalg import norm, solve

from scipy.linalg import solve_banded

import matplotlib.pyplot as plt




def updateTemperature(Er,
                      Cl,
                      Cg,
                      idx,
                      Tt,
                      Tb,
                      kappa = 1.0,
                      thcond = 'uniform'):

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





       
    # Constant vector

    nn = len(Er)
    
    B = np.ones( nn )

    B.fill(0)

    B[0] = Tb

    B[-1] = Tt
    



    # Matrix construction
    
    A = np.zeros( (nn, nn) )
    

    for i in range(nn):


        
        # First or last node
        
        if (i == 0) or (i == nn-1):

            A[i,i] = 1.


            
        # Inner nodes
            
        else:

            # Thermal conductivities
        
            lm = kappa

            li = kappa

            lp = kappa


            if i < idx:

                lm = lambdaRho( Cl[i-1], kappa, thcond )

                li = lambdaRho( Cl[i], kappa, thcond )

                lp = lambdaRho( Cl[i+1], kappa, thcond )

            elif i == idx:
            
                lm = lambdaRho( Cl[i-1], kappa, thcond )

                li = lambdaRho( 0.5*Cl[i] + 0.5*Cg[i], kappa, thcond )

                lp = lambdaRho( Cg[i+1], kappa, thcond )

            else:

                lm = lambdaRho( Cg[i-1], kappa, thcond )

                li = lambdaRho( Cg[i], kappa, thcond )

                lp = lambdaRho( Cg[i+1], kappa, thcond )


                
            # First/second order coefficients

            A[i,i-1] = lm
            
            A[i,i] = -lm - lp
            
            A[i,i+1] = lp
        

            

    # # Solucion del sistema    

    # Tr = solve(A,B)



    # Solucion usando banded

    ab = np.zeros((3,nn))

    for j in range(1,nn):    

        ab[0,j] = A[j-1,j]

        
    for j in range(nn):    

        ab[1,j] = A[j,j]


    for j in range(0,nn-1):    

        ab[2,j] = A[j+1,j]



    Tr = solve_banded((1,1),ab,B)


    
                
    return Tr
