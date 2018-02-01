import numpy as np


def analiticInterphase( er, cr ):

    xcross = 0.
    
    for i in range( len(er) - 1 ):
        
        if er[i] == er[i+1]:

            xcross = er[i]


    return xcross

        
