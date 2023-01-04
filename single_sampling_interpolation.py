#!/bin/env python3

import subprocess

import numpy as np
from pyscript import *

method = 'newton'
name = 'ethane'
trajectory = int(input(f'which trajectory should be calculated?'))
threshold_DCP_guess = 1e-1
ProcessMaxima = 'no'
max_steps = 100
make_extended_calculation = True

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
        with open(f'./trajectory-{trajectory}/sampling/fort.100') as reffile:
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
    with open(f'trajectory-{trajectory}/sampling/none.ami', 'w') as printfile:
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

    #creates .yml file
    with open(f'trajectory-{trajectory}/sampling/cluster.yml', 'w') as printfile:
        printfile.write(f'''--- # MaximaProcessing
MaximaProcessing:
  binaryFileBasename: none
  calculateSpinCorrelations: false
  shuffleMaxima: false
...''')
    
    #executes amolqc and Processmaxima
    with cd(f'trajectory-{trajectory}/sampling'):
        success = True
        try:
            run(f'amolqc none.ami')
        except subprocess.CalledProcessError:
            success = False
        if ProcessMaxima == 'yes':
            if success:
                run('ProcessMaxima cluster.yml')
            else:
                print(f'method none did not work for trajectory-{trajectory}')
    
    
def reading_coordinates(trajectory, calculation_type):
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
    
    
#interpolation starts here
n_elecs = reading_n_elecs()

mkdir(f'trajectory-{trajectory}/sampling')
cp(f'{name}.wf', f'trajectory-{trajectory}/sampling')

start_pos1 = DCP_guess_point(basin_minimum(trajectory, 'basin1'), basin_minimum(trajectory, 'basin2'), basin_outer_point(trajectory, 'basin1'), threshold_DCP_guess)
start_pos2 =  DCP_guess_point(basin_minimum(trajectory, 'basin1'), basin_minimum(trajectory, 'basin2'), basin_outer_point(trajectory, 'basin2'), threshold_DCP_guess)

a = start_pos1
single_point_none(trajectory, n_elecs, a)
start_pot1 = phi_value(trajectory)

b = start_pos2
single_point_none(trajectory, n_elecs, b)
start_pot2 = phi_value(trajectory)

list_potentials = []
list_position = []

for step in range(max_steps+1):
    intervall = 1/(max_steps)
    single_point_none(trajectory, n_elecs, (intervall * step * a + (1- intervall * step) * b ))
    pot = np.round(phi_value(trajectory), 5)
    position = np.round(intervall * step, 5)
    list_potentials.append(pot)
    list_position.append(position)
    
    #calculates distance from guess point to found DCP
    distance_to_DCP = np.linalg.norm(reading_coordinates(trajectory, method) - (intervall * step * a + (1- intervall * step) * b ))
    list_distance_to_DCP.append(distance_to_DCP)
    
#calculate distance of both positions
length_outer_basin_points = np.linalg.norm(start_pos1 - start_pos2)


    
    
with open(f'trajectory-{trajectory}/result/sampling.out', 'w') as printfile:
    for i_item, item in enumerate(list_potentials):
        printfile.write(f'{list_position[i_item]}a: {item} \n')
    printfile.write(f'maximum value at {list_position[list_potentials.index(max(list_potentials))]}: {max(list_potentials)}\n')
    printfile.write(f'list_position = {list_position}\n')
    printfile.write(f'list_potentials = {list_potentials}\n')
    printfile.write(f'list_distance_to_DCP = {list_distance_to_DCP}\n')
    printfile.write(f'norm = {length_outer_basin_points}\n')
print('sampling.out file was generated in folder result')
