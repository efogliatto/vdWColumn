import numpy as np

from numpy.linalg import norm, solve




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
        
            lm = 1.

            li = 1.

            lp = 1.


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

                



            # Second node

            if i == 1:

                A[i,i-1] = 3.*lm

                A[i,i]   = -4.*lm - lp

                A[i,i+1] = lm

                A[i,i+2] = lp
            

            # Last node

            elif i == nn-2:

                A[i,i-2] = lm
                
                A[i,i-1] = -3.*lp

                A[i,i]   = 4.*lp - lm

                A[i,i+1] = -lp


            # Inner nodes

            else:

                A[i,i-2] = lm
                                                    
                A[i,i]   = -lm - lp                

                A[i,i+2] = lp




    # Solucion del sistema    

    Tr = solve(A,B)


                
    return Tr
