import numpy as np

import matplotlib.pyplot as plt

import glob


# Find max position and half maximum limits using second derivative

def ddlimits( first, second, width = 0.1 ):

    
    # Discrete maximum (abs) position

    dpos = np.argmax( abs(first) )

    

    # find left limit
    
    aa = dpos-1

    while( (second[aa] <= second[aa+1])  and  (aa > 0) ):

        aa = aa - 1


    a = aa

    while(  ( second[a] <= second[aa]*width )   and   (a > 0)  ):

        a = a - 1    



        
    # find right limit
    
    bb = dpos+1

    while(   ( second[bb] >= second[bb-1] )  and   (bb < len(second)-1)   ):

        bb = bb + 1


    b = bb
        
    while(   ( second[b] >= second[bb]*width )   and   (b < len(second)-1)   ):

        b = b + 1    

    
    return dpos, a, b        
