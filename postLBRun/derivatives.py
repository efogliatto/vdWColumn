import numpy as np



# Compute derivatives of 1D array
    
def derivatives( field ):

    
    # Compute first and seconde derivatives

    dy  = np.zeros( len(field) )

    ddy = np.zeros( len(field) )    

    
    for i in range(1, len(dy) - 1):

        dy[i]  = field[i+1] - field[i-1]

        ddy[i] = 0.5  *  ( field[i+1] + field[i-1] - 2 * field[i] )


    return dy, ddy
