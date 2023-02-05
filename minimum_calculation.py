#!/bin/env python3

import subprocess

from pyscript import *
import re
import numpy as np


def reading_n_elecs():
    with open(f'trajectory-1-max.ref', 'r') as reffile:
        for line in reffile:
            if 'MAX:' in line:  
                line = reffile.readline()
                words = line.split()
                n_elecs = int(words[0])
    return n_elecs


def read_ref_file(reffile, coordinate_position):
    temp = reffile.upper()
    with open(f'trajectory-{trajectory}-{reffile}.ref', 'r') as reffile:
        R = []
        for line in reffile:
        if f'{coordinate_position} F({temp}):' in line:
            line = reffile.readline()
            for _ in range(n_elecs):
                line = reffile.readline()
                words = line.split()
                for word in words:
                    R.append(float(word))
            return np.array(R)
        
        
def read_coordinates(trajectory, calculation_type):
    if calculation_type == method:
        path = f'trajectory-{trajectory}/DCP_{method}/fort.100'
    elif calculation_type == 'stedes_eigvec':
        path = f'eigenvector_check/fort.100'
    with open(path) as reffile:
            R = []
            for line in reffile:
                if 'after minimize:' in line:
                    line = reffile.readline()
                    line = reffile.readline()
                    for _ in range(n_elecs):
                        line = reffile.readline()
                        words = line.split()
                        for word in words[1:]:
                            R.append(float(word))
    return(np.array(R))


def read_trajectory_ami(search):
    with open(f'trajectory.ami', 'r') as reffile:
        if search == 'count':
            for line in reffile:
               if 'count' in line:
                   count = int(re.search(r'\d+', line).group())
                   break
            return(count)
        if search == 'file':
            for line in reffile:
                if 'file' in line:
                    wavefunction = re.search(r'file=([\'"]?)(.+?)\.wf\1', line).group(2)
                    break
            return(wavefunction)


#script starts here
#definitions
count = read_trajectory_ami('count')
wavefunction = read_trajectory_ami('file')
n_elecs = read_n_elecs()

for trajectory in range(1, count+1):
    #checks if a second minimum was identified and if the second potetial is infty
    with open(f'trajectory-{trajectory}-max.ref') as reffile:
        found = False
        infty = False
        for line in reffile:
            if '2 F(MAX):' in line:
                found = True
                words = line.split()
                if words[4] == 'Infinity':
                    infty = True

    #calculates the minimum
    if found and not infty:
        #creates folderstructure
        mkdir(f'trajectory-{trajectory}')
        mkdir(f'trajectory-{trajectory}/result')
        cp(f'trajectory-start/result/cluster_start-out.yml', f'trajectory-{trajectory}/result')
        mkdir(f'trajectory-{trajectory}/minimum')
        cp(f'{wavefunction}.wf', f'trajectory-{trajectory}/minimum')

        #reads out the coordinates of the ending minimum
        R_max = read_ref_file('max', 2)

        #creates .ami file for Amolqc to make a single point run
        with open(f'trajectory-{trajectory}/minimum/stedes.ami', 'w') as sted_file:
            sted_file.write(f'''! seed for random number generation, not important
$gen(seed=101)
! reading the wave function file
$wf(read,file='{wavefunction}.wf')
$init_rawdata_generation()
$init_max_search(
step_size=0.1,
correction_mode=cut,
singularity_threshold=0.001,
method=steepest_descent,
verbose=2,
negative_eigenvalues=0,
eigenvalue_threshold=1e-10)
! setting the initial position
$init_walker(
free
''')
            for s in range(n_elecs):
                for t in R_max[s*3:s*3+3]:
                    sted_file.write(f'{t} ')
                sted_file.write('\n')
            sted_file.write(''')
$sample(create, size=1, single_point)
! maximize the walker
$maximize_sample()''')
       
        #creates .yml file
        with open(f'trajectory-{trajectory}/minimum/cluster.yml', 'w') as printfile:
            printfile.write(f'''--- # MaximaProcessing
MaximaProcessing:
  binaryFileBasename: stedes
  calculateSpinCorrelations: false
  shuffleMaxima: false
...''')

        #makes the amolqc run for the minimum single point calculation
        with cd(f'trajectory-{trajectory}/minimum'):
            run('amolqc stedes.ami')
            run('ProcessMaxima cluster.yml')
            cp('cluster-out.yml', '../result/cluster_end-out.yml')

        print(f'Minimum calculated for trajectory-{trajectory}')
        
    else: print(f'Minimum NOT calculated for trajectory-{trajectory}')   
