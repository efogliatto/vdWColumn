import numpy as np

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
