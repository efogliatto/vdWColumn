from .interphaseDensities import interphaseDensities


def latentHeat( T, a, b, R = 1 ):


    Tc = 8 * a / (27 * R * b)
    
    cl, cg = interphaseDensities( T / Tc )

    rl = cl / (3 * b)

    rg = cg / (3 * b)


    hl = (3./2. + 1./(1. - b*rl)) * R * T - 2 * a * rl

    hg = (3./2. + 1./(1. - b*rg)) * R * T - 2 * a * rg

    

    return hl, hg
