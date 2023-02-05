#!/bin/env python3

import subprocess

from pyscript import *
import re
import numpy as np


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
                    name = re.search(r'file=([\'"]?)(.+?)\.wf\1', line).group(2)
                    break
            return(name)
        
        
def read_n_elecs():
    with open(f'trajectory-1-max.ref', 'r') as reffile:
        line = reffile.readline()
        while 'MAX:' not in line:
            line = reffile.readline()
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
            

def interpolation(R_OBP1, R_OBP2):
    return( R_OBP1 + R_OBP2) / 2

    
def initial_guess_point(R_max1, R_max2, R_int, threshold_molecule):
    for l in range(n_elecs):
        norm = np.linalg.norm(R_max1[l*3:l*3+3]-R_max2[l*3:l*3+3])
        if norm < threshold_molecule:
            for r in R_max2[l*3:l*3+3]:
                printfile.write(f'{r} ')
        else:
            for t in R_int[l*3:l*3+3]:
                printfile.write(f'{t} ')
        printfile.write('\n')
        

def DCP_coordinates(trajectory, method):
    with open(f'trajectory-{trajectory}/DCP_{method}/fort.100') as reffile:
            R_DCP = []
            for line in reffile:
                if 'after minimize:' in line:
                    line = reffile.readline()
                    line = reffile.readline()
                    for _ in range(n_elecs):
                        line = reffile.readline()
                        words = line.split()
                        for word in words[1:]:
                            R_DCP.append(float(word))
    return(np.array(R_DCP))
        
    
def read_saddlepoint_calculation_in(search):
    with open('saddlepoint_calculation.in', 'r') as reffile:
        found = False
        for line in reffile:
            if search in line:
                found = True
                words = line.split()
                try:
                    float(words[1])
                    is_float = True
                except ValueError:
                    is_float = False
                if is_float:
                    return(float(words[1]))
                else:
                    return(words[1])
        #default values if not defined in .in file
        if not found:
            if search == 'threshold_DCP_guess':
                return(1e-1)
            if search == 'method':
                return('gradient_norm')
            if search  == 'deflection_factor':
                return(3e-3)
            

#script starts here
#reads out and defines count and n_elecs
n_elecs = read_n_elecs()
count = read_trajectory_ami('count')
name = read_trajectory_ami('file')
            
#reads saddlepoint_calculation.in file 
threshold_molecule = read_saddlepoint_calculation_in('threshold_DCP_guess')
method = read_saddlepoint_calculation_in('method')

#loop for trajectories starts here        
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

    if found and not infty:
    #calculates saddlepoints with method from input
        #makes folders
        mkdir(f'trajectory-{trajectory}/DCP_{method}')
        cp(f'{name}.wf', f'trajectory-{trajectory}/DCP_{method}')


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
singularity_threshold=0.0001,
''')
            #prints inputs for the chosen method
            if method == 'newton':
                printfile.write('method=newton,')
            elif method == 'none':
                printfile.write('method=none,')
            elif method == 'gradient_norm':
                printfile.write('convergence_value=5.0e-4,\n')
                printfile.write('method=bfgs,\n') 
                printfile.write('minimize_grad_norm,\n') 
                printfile.write('switch_step=50,')
            printfile.write('''
verbose=2,
negative_eigenvalues=-1,
eigenvalue_threshold=1e-10)
! setting the initial position
$init_walker(
free
''')
            #prints the initial guess point
            initial_guess_point(read_ref_file('max', 1), read_ref_file('max', 2),  interpolation(read_ref_file('traj', 1), read_ref_file('traj', 2)), threshold_molecule)
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

                
#makes a single point newton calculation, to obtain the order of the saddle point and the eigenvectors if gradient norm is the used method
        if method == 'gradient_norm':
            mkdir(f'trajectory-{trajectory}/DCP_{method}/newton_singlepoint')
            cp(f'{name}.wf', f'trajectory-{trajectory}/DCP_{method}/newton_singlepoint')
            
            with open(f'trajectory-{trajectory}/DCP_{method}/newton_singlepoint/newton_singlepoint.ami', 'w') as printfile:
                printfile.write(f'''! seed for random number generation, not important
$gen(seed=101)
! reading the wave function file
$wf(read,file='{name}.wf')
$init_rawdata_generation()
$init_max_search(
step_size=0.1,
correction_mode=none,
singularity_threshold=0.0001,
method=newton,
verbose=2,
negative_eigenvalues=-1,
eigenvalue_threshold=1e-10)
! setting the initial position
$init_walker(
free
''')
                R_DCP = DCP_coordinates(trajectory, method)
                for l in range(n_elecs):    
                    for t in R_DCP[l*3:l*3+3]:
                        printfile.write(f'{t} ')
                    printfile.write('\n')
                printfile.write(''')
$sample(create, size=1, single_point)
! maximize the walker
$maximize_sample()''')
            with cd(f'trajectory-{trajectory}/DCP_{method}/newton_singlepoint'):
                success = True
                try:
                    run('amolqc newton_singlepoint.ami')
                except subprocess.CalledProcessError:
                    success = False
                    
            if success:
                print(f'trajectory = {trajectory} newton_singlepoint was calculated')
            else:
                print(f'trajectory = {trajectory} newton_singlepoint could not be calculated')

    else:
        print(f'trajectory-{trajectory} could not be calculated')                                                     
