#!/bin/env python3

import subprocess

import numpy as np
from pyscript import *

method = 'newton'
name = 'ethane'
trajectory = int(input('which trajectory should be calculated?'))
threshold_DCP_guess = 1e-1

def reading_n_elecs():
    with open(f'trajectory-1-max.ref', 'r') as reffile:
        line = reffile.readline()
        while 'MAX:' not in line:
            line = reffile.readline()
        line = reffile.readline()
        words = line.split()
        n_elecs = int(words[0])
    return n_elecs


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


def basin_outer_point(trajectory, basin):
    with open(f'./trajectory-{trajectory}-traj.ref', 'r') as reffile:
        if basin == 'basin1':
            x = 1 
        elif basin == 'basin2':
            x = 2 
        R = []
        line = reffile.readline()
        for line in reffile: 
            if f'{x} F(TRAJ):' in line:
                line = reffile.readline()
                for _ in range(n_elecs):
                    line = reffile.readline()
                    words = line.split()
                    for word in words:
                        R.append(float(word))
    return np.array(R)


def basin_minimum(trajectory, basin):
    with open(f'./trajectory-{trajectory}-max.ref', 'r') as reffile:
        if basin == 'basin1':
            x = 1 
        elif basin == 'basin2':
            x = 2 
        R = []
        line = reffile.readline()
        for line in reffile: 
            if f'{x} F(MAX):' in line:
                line = reffile.readline()
                for _ in range(n_elecs):
                    line = reffile.readline()
                    words = line.split()
                    for word in words:
                        R.append(float(word))
    return np.array(R)


def DCP_guess_point(R_max1, R_max2, basin_outer_point, threshold_DCP_guess):
    R = []
    for l in range(n_elecs):
        norm = np.linalg.norm(R_max1[l*3:l*3+3]-R_max2[l*3:l*3+3])
        if norm < threshold_DCP_guess:
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
    

#interpolation starts here
n_elecs = reading_n_elecs()

mkdir(f'trajectory-{trajectory}/temp')
cp(f'{name}.wf', f'trajectory-{trajectory}/temp')

start_pot1 =  DCP_guess_point(basin_minimum(trajectory, 'basin1'), basin_minimum(trajectory, 'basin2'), basin_outer_point(trajectory, 'basin1'), threshold_DCP_guess)
start_pot2 =  DCP_guess_point(basin_minimum(trajectory, 'basin1'), basin_minimum(trajectory, 'basin2'), basin_outer_point(trajectory, 'basin2'), threshold_DCP_guess)

a = start_pot1
single_point_none(trajectory, n_elecs, a)
print(phi_value(trajectory))

b = start_pot2
single_point_none(trajectory, n_elecs, b)
print(phi_value(trajectory))


def f(x, trajectory):
    single_point_none(trajectory, n_elecs, x)
    return(phi_value(trajectory))


def gss(f, a, b, tol=1e-1):

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

gr = (np.sqrt(5) + 1) / 2 

R_new = gss(f, a, b)

print(R_new)
print(f(R_new, trajectory))
                           
