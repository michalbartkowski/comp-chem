# 2018-05-19 Michal B
'''
This script batch converts molecular coordinate (*.xyz) files
into gaussian09 input (*.com) files. Specify *.com settings
under the INPUTS section of this script '''

import numpy as np
import os
import sys
import glob

# INPUTS ----------
theory = '#n B3LYP/6-31G(d) SP'
charge = 0
multiplicity = 1
# -----------------


# Extracts data from *.xyz file
def get_data(file):


    with open(file, "r") as fo:
        lines = fo.readlines()
    data = [i.split() for i in lines[2:]]
    #data = np.array(data)
    #data.insert(0,theory)
    return data

if __name__ == '__main__':	
    
    # Work on every .xyz file in the batch folder specified
    os.chdir("xyz/")
    if not os.path.exists("gaussian_input"):
        os.makedirs("gaussian_input")
    for file in glob.glob("*.xyz"):
        # Removes path and extension; returning just the name
        filename = os.path.splitext(os.path.basename(file))[0]

        data = get_data(file)

        data.insert(0, [theory,'','',''])
        data.insert(1, ['','','',''])
        data.insert(2, [' ' + str(filename),'','',''])
        data.insert(3, ['','','',''])
        data.insert(4, [charge, multiplicity,'',''])

        
        #----WRITE DATA OUTPUT----
        np.savetxt('gaussian_input/' + str(filename) + '.com', data, fmt="%s")