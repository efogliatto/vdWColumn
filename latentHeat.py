from .interphaseDensities import interphaseDensities


def latentHeat( T, a, b, R = 1 ):


    Tc = 8 * a / (27 * R * b)
    
    cl, cg = interphaseDensities( T / Tc )

    rl = 3 * b * cl

    rg = 3 * b * cg


    hl = (3./2. + 1./(1. - b*rl)) * R * T - 2 * a * rl

    hg = (3./2. + 1./(1. - b*rg)) * R * T - 2 * a * rg


    return hl, hg
