import numpy as np

import matplotlib.pyplot as plt

from .derivatives import derivatives

from .fwhm import fwhm

from .ddlimits import ddlimits



def interphase( field, width = 0.5, fw = False ):


    # Compute derivatives

    dy,ddy = derivatives( field )

    

    dpos, ll, rl = 0, 0, 0

    # Discrete maximum (abs) derivative position
    
    if fw == True:

        dpos, ll, rl = fwhm( abs(dy), width )

        
    else:

        dpos, ll, rl = ddlimits( dy, ddy )
        



    xcross = 0

    if (   (ll > 0)   and   (rl < len(field))   ):


        # Slice second derivative and get zero crossing

        dcross = 0

        for i in range(dpos-5, dpos+5):       
        
            if ddy[i] * ddy[i+1] <= 0.:

                dcross = i



        # Get continuous position

        a = ddy[dcross+1] - ddy[dcross]

        b = ddy[dcross] - a * dcross


        xcross = - b / a
    
    

        

    return xcross, ll, rl
