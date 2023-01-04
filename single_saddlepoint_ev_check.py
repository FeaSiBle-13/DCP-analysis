#!/bin/env python3

import subprocess

from pyscript import *
import numpy as np

deflection_factor = 1e-3
method = 'newton'
trajectory = int(input('which trajectory do you mean?'))
name = 'ethane'
threshold = 1e-2

def reading_n_elecs():
    with open(f'trajectory-1-max.ref', 'r') as reffile:
        line = reffile.readline()
        while 'MAX:' not in line:
            line = reffile.readline()
        line = reffile.readline()
        words = line.split()
        n_elecs = int(words[0])
    return n_elecs


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


def minimum(trajectory, x):
    with open(f'trajectory-{trajectory}-max.ref', 'r') as reffile:
        R = []
        line = reffile.readline()
        while f'{x} F(MAX):' not in line:
            line = reffile.readline()
        line = reffile.readline()
        for _ in range(n_elecs):
            line = reffile.readline()
            words = line.split()
            for word in words:
                R.append(float(word))
        return np.array(R)


def read_eigenvector(trajectory):
    with open(f'trajectory-{trajectory}/DCP_newton/fort.100') as reffile:
        R = []
        for line in reffile:
            found = True
            if 'hessian eigenvalues and -vectors' in line:
                line = reffile.readline()
                for _ in range(n_elecs):
                    line = reffile.readline()
                    coordinates = line.split()
                    for coordinate in coordinates[1:]:
                        R.append(float(coordinate))
    return np.array(R)


def deflection_saddlepoint(eigenvector, saddlepoint, deflection_factor):
    return(eigenvector * deflection_factor + saddlepoint)


def stepest_descent(trajectory, molecule, R):
    cp('ethane.wf', f'eigenvector_check')
    with open(f'eigenvector_check/stedes.ami', 'w') as printfile:
        printfile.write(f'''! seed for random number generation, not important
$gen(seed=101)
! reading the wave function file
$wf(read,file='{name}.wf')
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
            for t in R[s*3:s*3+3]:
                printfile.write(f'{t} ')
            printfile.write('\n')
        printfile.write(''')
$sample(create, size=1, single_point)
! maximize the walker
$maximize_sample()''')


    with open(f'eigenvector_check/cluster.yml', 'w') as printfile:
        printfile.write(f'''--- # MaximaProcessing
MaximaProcessing:
  binaryFileBasename: stedes
  calculateSpinCorrelations: false
  shuffleMaxima: false
...''')

    #makes the amolqc run for the minimum single point calculation
    with cd(f'eigenvector_check'):
        run('amolqc stedes.ami')
    return(print(f'minimum for eigenvector_check for trajectory-{trajectory} was calculated'))


def compare_position(R1, R2, threshold):
    norm = np.linalg.norm(R1 - R2)
    if norm <= threshold:
        return(True)
    else:
        return(False)


#program starts here
n_elecs = reading_n_elecs()

mkdir(f'eigenvector_check')

eigenvec = read_eigenvector(trajectory)
saddlepoint = reading_coordinates(trajectory, method)


#compares if stepest descent minimization ends in the starting or ending minimum
list_minimized_deflection = []

for m in range(1, 3):
    stepest_descent(trajectory, name, deflection_saddlepoint(eigenvec, saddlepoint, (-1) ** m * deflection_factor))
    list_minimized_deflection.append(reading_coordinates(trajectory, 'stedes_eigvec'))
    for obj in ls(f'eigenvector_check'):
        rm(f'eigenvector_check/{obj}')


list_found_min = [False, False]
for m in range(1, 3):
    for vectors in list_minimized_deflection:
        same = compare_position(minimum(trajectory, m), vectors, threshold)
        if same:
            list_found_min[m-1] = True
            break

if list_found_min[0] and list_found_min[1]:
    print('the DCP lies between the starting minimum and the second minimum of the trajectory')
else:
    print('the DCP lies NOT between the starting minimum and the second minimum of the trajectory')

