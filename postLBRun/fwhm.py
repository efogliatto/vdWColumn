import numpy as np



# Find max position and full widht half maximum limits. Specially for density derivative

def fwhm( field, width = 0.5 ):

    
    # Discrete maximum (abs) position

    dpos = np.argmax( abs(field) )


    # find left limit
    
    a = dpos

    while( field[a] >= field[dpos]*width ):

        a = a - 1



    # find right limit
    
    b = dpos

    while( field[b] >= field[dpos]*width ):

        b = b + 1



    return dpos, a, b
