#!/bin/env python3

import subprocess

from pyscript import *
import re
import numpy as np


def reading_n_elecs():
    with open(f'trajectory-1-max.ref', 'r') as reffile:
        line = reffile.readline()
        while 'MAX:' not in line:
            line = reffile.readline()
        line = reffile.readline()
        words = line.split()
        n_elecs = int(words[0])
    return n_elecs


def starting_maximum(trajectory):
    with open(f'trajectory-{trajectory}-max.ref', 'r') as reffile:
        R_max = []
        line = reffile.readline()
        while '1 F(MAX):' not in line:
            line = reffile.readline()
        line = reffile.readline()
        for _ in range(n_elecs):
            line = reffile.readline()
            words = line.split()
            for word in words:
                R_max.append(float(word))
        return np.array(R_max)


def ending_maximum(trajectory):
    with open(f'trajectory-{trajectory}-max.ref', 'r') as reffile:
        R_max = []
        line = reffile.readline()
        while '2 F(MAX):' not in line:
            line = reffile.readline()
        line = reffile.readline()
        for _ in range(n_elecs):
            line = reffile.readline()
            words = line.split()
            for word in words:
                R_max.append(float(word))
        return np.array(R_max)
        print(R_max)


def basin_left(trajectory):
    with open(f'./trajectory-{trajectory}-traj.ref', 'r') as reffile:
        R_i = []
        line = reffile.readline()
        while '1 F(TRAJ):' not in line:
            line = reffile.readline()
        line = reffile.readline()
        for _ in range(n_elecs):
            line = reffile.readline()
            words = line.split()
            for word in words:
                R_i.append(float(word))
    return np.array(R_i)


def basin_enter(trajectory):
    with open(f'./trajectory-{trajectory}-traj.ref', 'r') as reffile:
        R_j = []
        line = reffile.readline()
        while '2 F(TRAJ):' not in line:
            line = reffile.readline()
        line = reffile.readline()
        for _ in range(n_elecs):
            line = reffile.readline()
            words = line.split()
            for word in words:
                R_j.append(float(word))
    return np.array(R_j)
 
def DCP_guess_point(R_max1, R_max2, threshold_DCP_guess):
    for l in range(n_elecs):
        norm = np.linalg.norm(R_max1[l*3:l*3+3]-R_max2[l*3:l*3+3])
        if norm < threshold_DCP_guess:
            for r in R_max2[l*3:l*3+3]:
                printfile.write(f'{r} ')
        else:
            for t in R_int[l*3:l*3+3]:
                printfile.write(f'{t} ')
        printfile.write('\n')


#reads count from .ami
with open(f'trajectory.ami', 'r') as ami_file:
    notdone_count = True
    notdone_file = True
    for line in ami_file:
        if 'count' in line and notdone_count:
            count = int(re.search(r'\d+', line).group())
            notdone_count = False
        elif 'file' in line and notdone_file:
            name = re.search(r'file=([\'"]?)(.+?)\.wf\1', line).group(2)
            notdone_file = False
            
#reads n_elecs
n_elecs = reading_n_elecs()

#reads saddlepoint_calculation.in file (would be nicer here wirth regular expressions)
    with open('saddlepoint_calculation.in', 'r') as reffile:
        for line in reffile:
            if 'threshold_DCP_guess' in line:
                words = line.split()
                threshold_DCP_guess = int(words[1])
            else:
                threshold_DCP_guess = 1e-1
                
            if 'method' in line:
                words = line.split()
                method = words[1]
            else: 
                method = 'gradient_norm'
        
#loop for trajectories starts here        
for trajectory in range(1, count+1):
    with open(f'trajectory-{trajectory}-max.ref') as reffile:
        found = False
        infty = False
        for line in reffile:
            if '2 F(MAX):' in line:
                found = True
                words = line.split()
                if words[4] == 'Infinity':
                    infty = True

    if found and not infty:
        
    #calculates saddlepoints with method input
        #makes folders
        mkdir(f'trajectory-{trajectory}/DCP_{method}')
        for obj in ls(f'trajectory-source/DCP'):
            cp(f'trajectory-source/DCP/{obj}', f'trajectory-{trajectory}/DCP_{method}')


        #interpolation of both -traj- vectors
        R_int = ( basin_left(trajectory) + basin_enter(trajectory) ) / 2


        #creates the {method}.ami with the interpolated vectors and the input with the fixed e positions
        with open(f'trajectory-{trajectory}/DCP_{method}/{method}.ami', 'w') as printfile:
            printfile.write(f'''! seed for random number generation, not important
$gen(seed=101)
! reading the wave function file
$wf(read,file='{name}.wf')
$init_rawdata_generation()
$init_max_search(
step_size=0.1,
correction_mode=none,
singularity_threshold=0.0001,''')
             #continue here
            if method == 'newton':
                printfile.write('method=newton')
            elif method == 'none':
                printfile.write('method=none')
            elif method == 'gradient_norm':
                printfile.write('minimize_gradient_norm,\n') 
                printfile.write('switch_step=50,\n')
            printfile.write('\n')
            printfile.write('''
verbose=2,
negative_eigenvalues=-1,
eigenvalue_threshold=1e-10)
! setting the initial position
$init_walker(
free
''')
            DCP_guess_point(starting_maximum(trajectory), ending_maximum(trajectory), threshold_DCP_guess)
            printfile.write(''')
$sample(create, size=1, single_point)
! maximize the walker
$maximize_sample()''')

#creates .yml file
        with open(f'trajectory-{trajectory}/DCP_{method}/cluster.yml', 'w') as printfile:
            printfile.write(f'''--- # MaximaProcessing
MaximaProcessing:
  binaryFileBasename: {method}
  calculateSpinCorrelations: false
  shuffleMaxima: false
...''')
                             
#executes amolqc and Processmaxima
        with cd(f'trajectory-{trajectory}/DCP_{method}'):
            success = True
            try:
                run(f'amolqc {method}.ami')
            except subprocess.CalledProcessError:
                success = False

            if success:
                run('ProcessMaxima cluster.yml')
                cp('cluster-out.yml', f'../result/cluster_{method}_trans-out.yml')
            else:
                print(f'method {method} did not work for trajectory-{trajectory}')


        #executesEVanalysis
        if success:
            with cd(f'trajectory-{trajectory}/result'):
                try:
                    run(f'/home/theochem/Scripts/EVanalysis/EVanalysis.py ../minimum/{name}.wf cluster_{method}_trans-out.yml')
                except subprocess.CalledProcessError:
                    pass


    else:
        print(f'trajectory-{trajectory} could not be calculated')                                                     


