import math


# Reduced saturation pressure

def P_sat_r ( cg, cl ):

    """Saturation pressure
       
       Arguments
       cg: reduced concentration of vapor phase
       cl: reduced concentration of liquid phase"""
    
    return cg * cl * (3.0 - (cg + cl))
