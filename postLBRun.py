import numpy as np

import matplotlib.pyplot as plt

import glob




# Load 1D profile from scalar field

def scalarProfile( folderName, fieldName, step = 3, offset = 0, time = 'latest' ):

    if folderName[-1] == '/':

        fname = folderName + 'processor0/'

    else:
            
        fname = folderName + '/processor0/'


        
    # Check for latest time

    if time == 'latest':

            
        flist = glob.glob( fname + '*' )

        
        # Trim name

        tlist = []

        for a in flist:
            
            if a[len(fname):].isdigit():

                tlist.append(  int( a[len(fname):] )  )
        

        fname = fname + '{}'.format( sorted(tlist)[-1] ) + '/' + fieldName 



    else:

        fname = fname + '/' + time + '/' + fieldName


    
    

    fileField = np.loadtxt( fname, unpack = True)[1:]
    


    # Slice info

    idx = range( offset, len(fileField), step )
    
    field = np.zeros( len(idx) )
    
    for i,j in enumerate(idx):

        field[i] = fileField[j]



        

    return field
    




# Compute derivatives of 1D array
    
def derivatives( field ):

    
    # Compute first and seconde derivatives

    dy  = np.zeros( len(field) )

    ddy = np.zeros( len(field) )    

    
    for i in range(1, len(dy) - 1):

        dy[i]  = field[i+1] - field[i-1]

        ddy[i] = 0.5  *  ( field[i+1] + field[i-1] - 2 * field[i] )


    return dy, ddy






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






def analiticInterphase( er, cr ):

    xcross = 0.
    
    for i in range( len(er) - 1 ):
        
        if er[i] == er[i+1]:

            xcross = er[i]


    return xcross

        
