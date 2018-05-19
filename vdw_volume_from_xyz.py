''' 2018-05-04 -- Michal B

Script for calculating the vdW volume of n atoms by Monte Carlo method.
The script requires *.xyz files as input inside an \input directory, the
output is generated in a comma delimited format. The desired standard
deviation should be specified under the DESIRED_SD variable.
The vdw_rule dictionary should be updated with additional vdW radii if 
necessary.
'''

import multiprocessing as mp
import numpy as np
import time
import random
import os
import sys
import glob

# Determines whether the point x, y, z belongs to at least one of 
# the spheres. Calculates the dist squared of the point to the xyz 
# position of each atom, and checks if the point is within the vdw
# radii of said atom
def check_point(x,y,z,listatom):
    for i in range(listatom.size // 4):
        dist_sqr = (x - listatom[i,0]) ** 2 \
                 + (y - listatom[i,1]) ** 2 \
                 + (z - listatom[i,2]) ** 2
        if dist_sqr < ( listatom[i,3]) ** 2:
            return True

# Extracts data from *.xyz file
def get_data(file):
    #time_tt = time.time()
    vdw_values, xyz_values = [], []

    with open(file, "r") as fo:
        lines = fo.readlines()
    lines = [i.strip() for i in lines]
    # Dictionary for assigning vdw radii to respective atoms
    vdw_rule = {'H': 1.20, 'C' : 1.70, 'N' : 1.55, 'O': 1.52,
                'F': 1.47, 'Br': 1.85, 'Cl': 1.75, 'T': 1.00, 'Zn': 1.39}

    # Create array of vdw radii corresponding to their
    # respective atoms (using a dictionary: vdw_rule)
    for i in lines[2:]:
        if i.split()[0] in vdw_rule:
            vdw_values.append(vdw_rule[i.split()[0]])
    vdw_values = np.array(vdw_values).astype(float)

    # Create array of xyz coordinates
    #for i in lines[2:]:
    #   xyz_values.append(i.split()[1:])
    xyz_values = [i.split()[1:] for i in lines[2:]]
    xyz_values = np.array(xyz_values).astype(float)

    # merge coordinate array with respective vdw radii array
    listatom = np.insert(xyz_values, 3, vdw_values, axis=1)
    #print('Time to extract .xyz:', time.time()-time_tt)
    return listatom

# ----STANDARD DEVIATION----
def get_volume_sd(iterations, volumes, volume_mean):
    stdev = 0.0
    for i in range(0, iterations):
        stdev += (volumes[i] - volume_mean) ** 2
    stdev = (stdev / iterations) ** (0.5)
    return stdev

# ----AVERAGE VOLUME----
def get_volume_mean(iterations, volumes):
    volume_mean = 0.0
    for i in range(0, iterations):
        volume_mean = volume_mean + volumes[i]
    volume_mean /= iterations
    return volume_mean

def get_time(t):
    t_current = time.time() - t
    if  t_current <= 60:
        t_current = str(round(t_current, 1)) + ' s'
    else:
        t_current /= 60
        t_current = str(round(t_current, 1)) + ' min'	
    return t_current	

def main_process(listatom, points, xmin, xmax, ymin, ymax, zmin,
                 zmax, x, CPUs, output):
    random.seed()
    M = [0,0]
    for i in range(int(points/CPUs)):
        xx = random.uniform(xmin, xmax)
        xy = random.uniform(ymin, ymax)
        xz = random.uniform(zmin, zmax) 
        if check_point(xx, xy, xz, listatom): M[0] += 1
        else: M[1] += 1
    output.put(M)

if __name__ == '__main__':	
    
    DESIRED_SD = 3
    CPUs = 8
    iterations = 10

    # Work on every .xyz file in the batch folder specified
    os.chdir("input/")
    for file in glob.glob("*.xyz"):
        
        # Get xyz and vdw radii data
        listatom = get_data(file)

        # Declare array of 0's, which will be substituted for 
        # volume data from each calculation
        volumes = np.zeros(iterations)
        
        # Removes path and extension; returning just the name
        filename = os.path.splitext(os.path.basename(file))[0]

        # Declare cuboid vertices as floats
        minx, miny, minz =  999.0,  999.0,  999.0	
        maxx, maxy, maxz = -999.0, -999.0, -999.0	
    
        # Aligns cuboid constrainsts to xyz coordinates
        for nat in range(0, listatom.size // 4):
            if (listatom[nat, 0] < minx): minx = listatom[nat, 0]
            if (listatom[nat, 1] < miny): miny = listatom[nat, 1]
            if (listatom[nat, 2] < minz): minz = listatom[nat, 2]
            if (listatom[nat, 0] > maxx): maxx = listatom[nat, 0]
            if (listatom[nat, 1] > maxy): maxy = listatom[nat, 1]
            if (listatom[nat, 2] > maxz): maxz = listatom[nat, 2]                
    
        # Adds wiggle room, has to be larger than vdw radii
        minx -= 1.85
        miny -= 1.85
        minz -= 1.85
        maxx += 1.85
        maxy += 1.85
        maxz += 1.85       
    
        # Total cuboid volume
        v_cube = (maxx - minx) * (maxy - miny) * (maxz - minz)

        # Calculates random points required to achieve desired SD
        points = int(1107.3 * DESIRED_SD ** -1.954 * v_cube / \
                     iterations)
    
        # Print some statistics
        print('\n|  Data in the form: x, y, z, vdw radii:', 
              '\n', listatom,
              '\n|  Molecule:', filename,
              '\n|  Atoms:', listatom.size // 4,
              '\n|  Cuboid Volume:', round(v_cube, 3),
              '\n|  Total Points per Iteration:', points,
              '\n|  Iterations:', iterations,
              '\n|  CPUs Used:', CPUs,
              '\n|  Accuracy:', round(iterations * (points / \
                                      v_cube), 3))

        # Process to run for each iteration
        for v in range(0, iterations):
            # Define an output queue
            output = mp.Queue()  
    
            # Declare some timestamps
            if v == 0: time_init = time.time()
            time_it = time.time()
    
            # Setup a list of processes that we want to run
            processes = [mp.Process(target=main_process, 
                                    args=(listatom, points, minx, 
                                    maxx, miny, maxy, minz, maxz, 
                                    x, CPUs, output,)) 
                        for x in range(CPUs)]
            
            # Run processes
            for p in processes:
                p.start()
            
            # Exit the completed processes
            for p in processes:
                p.join() 
            	## print(pr, pr.is_alive())
    
            # Get process results from the output queue
            results = [output.get() for pr in processes]
            ## print(results)
            
            totin, totout = 0, 0
            for i in range(0, len(results)):
                totin = totin + results[i][0]
                totout = totout + results[i][1]
    
            # Volume calculated based on ratio of points in to out
            volumes[v] = v_cube * totin / (totin + totout)
    
            # Print statistics for each iteration
            print('\n---- ITERATION', v + 1, '- [' + \
                  get_time(time_it) + '] ----',
                  '\n  Volume:', round(volumes[v], 3),
                  '\n  Points In, Out:', totin, totout)
            
        # Print final results and statistics
        volume_mean = get_volume_mean(iterations, volumes)
        stdev = get_volume_sd(iterations, volumes, volume_mean)
        print('\n|  Mean Volume:', round(volume_mean, 3),
              '\n|  Population SD:', round(stdev, 3), 
              '\n|  Total Calc. Time:', get_time(time_init))	
        
        # ----WRITE DATA OUTPUT----
        #with open('output_SD_'+str(DESIRED_SD)+'.csv', 'a') as f:
        with open('output.csv', 'a') as f:
            # converts xyz/Broo.xyz to Broo
            f.write(str(os.path.splitext(os.path.basename(file))[0]) 
            + ',' + str(volume_mean) 
            + ',' + str(stdev)
            + ',' + str(v_cube)
            + ',' + str(points)
            + ',' + str(iterations * (points / v_cube))
            + ',' + str(round(time.time() - time_init,3)) 
            + ',' + str(iterations) + str('\n'))