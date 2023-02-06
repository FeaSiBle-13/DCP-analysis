#!/bin/env python3

import subprocess

import numpy as np
from pyscript import *


def read_n_elecs():
    with open(f'trajectory-1-max.ref', 'r') as reffile:
        for line in reffile:
            if 'MAX:' in line:  
                line = reffile.readline()
                words = line.split()
                n_elecs = int(words[0])
    return n_elecs


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


def read_DCP_analysis_in(search):
    with open('DCP_analysis.in', 'r') as reffile:
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
            if search == 'threshold_electrons':
                return(1e-3)
            if search == 'threshold_molecule':
                return(1e-1)
            if search == 'method':
                return('gradient_norm')
            if search  == 'deflection_factor':
                return(3e-3)
            if search  == 'compare_mode':
                return('electron_wise')



def phi_value(trajectory):
        with open(f'./trajectory-{trajectory}/temp/fort.100') as reffile:
            found = False
            for line in reffile:
                if 'Phi:' in line:
                    found = True
                    words = line.split()
                    phi = float(words[1])
                    break
            if not found:
                print('Phi from trajectory-start was not found')
        return(phi)


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


def initial_guess_point(R_max1, R_max2, basin_outer_point, threshold_DCP_guess):
    R = []
    for l in range(n_elecs):
        norm = np.linalg.norm(R_max1[l*3:l*3+3]-R_max2[l*3:l*3+3])
        if norm < threshold_electrons:
            for r in R_max2[l*3:l*3+3]:
                R.append(r)
        else:
            for t in basin_outer_point[l*3:l*3+3]:
                R.append(t)
    return(np.array(R))


def single_point_none(trajectory, n_elecs, R_elecs):
    with open(f'trajectory-{trajectory}/temp/vectors_potential.ami', 'w') as printfile:
        printfile.write(f'''! seed for random number generation, not important
$gen(seed=101)
! reading the wave function file
$wf(read,file='{name}.wf')
$init_rawdata_generation()
$init_max_search(
step_size=0.1,
correction_mode=none,
singularity_threshold=0.0001,
method=none,
verbose=2,
negative_eigenvalues=-1,
eigenvalue_threshold=1e-10)
! setting the initial position
$init_walker(
free
''')
        R = R_elecs
        for l in range(n_elecs):    
            for t in R[l*3:l*3+3]:
                printfile.write(f'{t} ')
            printfile.write('\n')
        printfile.write(''')
$sample(create, size=1, single_point)
! maximize the walker
$maximize_sample()''')
    with cd(f'trajectory-{trajectory}/temp'):
        success = True
        try:
            run(f'amolqc vectors_potential.ami')
        except subprocess.CalledProcessError:
            success = False
    if not success:
        print('Potential for golden section search could not be obtained')
        
        
def f(x, trajectory):
    single_point_none(trajectory, n_elecs, x)
    return(phi_value(trajectory))


def golden_sec_search(f, a, b, tol=1e-1):

    c = b - (b - a) / gr
    d = a + (b - a) / gr

    for step in range(0, 30):
        if abs(np.linalg.norm(b - a)) > tol:
            if f(c, trajectory) > f(d, trajectory):  # f(c) > f(d) to find the maximum
                b = d 
            else:
                a = c 
        else:
            break

        # We recompute both c and d here to avoid loss of precision which may lead to incorrect results or infinite loop
        c = b - (b - a) / gr
        d = a + (b - a) / gr
        print(f'step ={step+1}')
        print(f'a = {a}')
        print(f'b = {b}')
        print(f'c = {c}')
        print(f'd = {d}\n')
    return ((b + a) / 2)
    

#script starts here
#reads n_elecs from .max file and wavefunction from trajectory.ami
n_elecs = read_n_elecs()
name = read_trajectory_ami('file')

threshold_electrons = read_DCP_analysis_in('threshold_electrons')
method = read_DCP_analysis_in('method')

#manual input of the trajectory
trajectory = int(input('which trajectory should be calculated?'))

#creates file structure
mkdir(f'trajectory-{trajectory}/temp')
cp(f'{name}.wf', f'trajectory-{trajectory}/temp')

#makes the interpolation for the initial guess point
a =  initial_guess_point(read_ref_file('max', 1), read_ref_file('max', 2), read_ref_file('traj', 1), threshold_DCP_guess)
b =  initial_guess_point(read_ref_file('max', 1), read_ref_file('max', 2), read_ref_file('traj', 2), threshold_DCP_guess)

#makes single point none calculation to obtain the potential
single_point_none(trajectory, n_elecs, a)
print(phi_value(trajectory))

single_point_none(trajectory, n_elecs, b)
print(phi_value(trajectory))

#defines the golden ratio
gr = (np.sqrt(5) + 1) / 2 

#makes the golden section search
R_new = golden_sec_search(f, a, b)

print(R_new)
print(f(R_new, trajectory))
                           
