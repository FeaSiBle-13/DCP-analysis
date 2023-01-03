#!/bin/env python3

import numpy as np
from pyscript import *
import re
from sigfig import round

list_failure_DCP = []
list_failure_statistics = []
list_failure_trajectories = []
list_failure_phi_values = []
list_failure_order = []


def read_trajectory.ami(search):
    with open(f'trajectory.ami', 'r') as reffile:
        if search == 'count':
            for line in reffile:
               if 'count' in line and notdone_count:
                   count = int(re.search(r'\d+', line).group())
                   break
        if search == 'file':
            for line in reffile:
                if 'file' in line and notdone_file:
                name = re.search(r'file=([\'"]?)(.+?)\.wf\1', line).group(2)
                break
        
        
    return(count)


def phi_value(x, trajectory):
    if x == 0:
        path = f'trajectory-start/minimum/fort.100'
    elif x == 1:
        path = f'trajectory-{trajectory}/DCP_{method}/fort.100'
    elif x == 2:
        path = f'trajectory-{trajectory}/minimum/fort.100'
    with open(path) as reffile:
        found = False
        for line in reffile:
            if 'Phi:' in line:
               found = True
               words = line.split()
               phi = float(words[1])
        if not found:
            print(f'Phi from {path} was not found')
    return phi


def reading_order(trajectory):
    if method == 'gradient_norm':
        path = f'trajectory-{trajectory}/DCP_{method}/newton_singlepoint/fort.100'
    else:
        path = f'trajectory-{trajectory}/DCP_{method}/fort.100'
    with open(path) as reffile:    
        order = 0
        found = False
        for line in reffile:
            if 'hessian eigenvalues and -vectors:' in line:
                found = True
                for _ in range(10):
                    line = reffile.readline()
                    words = line.split()
                    if float(words[1]) < 0:
                        order += 1
                    else:
                        break
                    for _ in range(n_elecs):
                        line = reffile.readline()
        if not found:
            print(f'hessian eigenvalues and -vectors: was not found in file trajectory-{trajectory}/DCP_{method}/fort.100')
    return order


def reading_n_elecs():
    with open(f'trajectory-1-max.ref', 'r') as reffile:
        for line in reffile:
            if 'MAX:' in line:  
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


def compare_position(R1, R2, threshold_molecule):
    norm = np.linalg.norm(R1 - R2)
    if norm <= threshold_molecule:
        return(True)
    else:
        return(False)
    
    
def ev_deflection_check():
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
        return('saddlepoint_of_adjacent_minima')
    else:
        return('no_saddlepoint_of_adjacent_minima')
    

def compare_saddlepoints(R_new, trajectory):
    list_DCP = []
    list_phi_values = []
    list_enumerate_DCP =[]
    list_statistics = []
    list_trajectories = []
    list_order = []
    
    found = False
    for i_DCP, R_DCP in enumerate(list_DCP):
        norm = np.linalg.norm(R_new - R_DCP)
        if norm <= threshold_molecule:
            list_statistics[i_DCP] += 1
            list_trajectories[i_DCP] += f', {trajectory}'
            found = True
            break
    if not found:
        list_DCP.append(R_new)
        list_statistics.append(1)
        list_trajectories.append(f'{trajectory}')
        list_enumerate_DCP.append(len(list_enumerate_DCP)+1)
        list_phi_values.append(round(abs(phi_value(0, trajectory)-phi_value(1, trajectory)), sigfigs = 3 ))
        list_order.append(f'{reading_order(trajectory)}. order')
    list_compared[list_DCP, list_statistics, list_trajectories, list_enumerate_DCP, list_phi_values, list_order]
    return(list_compared)


def no_saddlepoint(reason):
    found = False
    for i_no_DCP, no_DCP in enumerate(list_failure_DCP):
        if no_DCP == reason:
            list_failure_statistics[i_no_DCP] += 1
            list_failure_trajectories[i_no_DCP] += f', {trajectory}'
            found = True
            break
    if not found:
        list_failure_DCP.append(reason)
        list_failure_statistics.append(1)
        list_failure_trajectories.append(f'{trajectory}')
        list_failure_phi_values.append('none')
        list_failure_order.append('none')
        

def runtime():
    with open('trajectory.amo', 'r') as reffile:
        for line in reffile:
            if 'wall clock time for run' in line:
                words = line.split()
                runtime = words[7]
    return runtime


def runtime_hours(runtime):
    time = runtime.split(':')
    runtime_hours = 0
    for i_item, item in enumerate(time):
        runtime_hours += float(item) / (60 ** float(i_item))
    return np.round(runtime_hours, decimals = 3)


#script starts here
#reads out and defines count and n_elecs
n_elecs = reading_n_elecs()
count = read_trajectory.ami('count')
name = read_trajectory.ami('file')
            
#reads saddlepoint_calculation.in file (would be nicer here with regular expressions)
with open('saddlepoint_calculation.in', 'r') as reffile:
    for line in reffile:
        if 'threshold_DCP_guess' in line:
            words = line.split()
            threshold_molecule = float(words[1])
        else:
            threshold_molecule = 1e-1
            
        if 'method' in line:
            words = line.split()
            method = words[1]
        else: 
            method = 'gradient_norm'            
            
            
#starts evaluation
mkdir('eigenvector_check')

for trajectory in range(1, count + 1):
    #checks if DCP was calculated
    with open(f'trajectory-{trajectory}-max.ref') as reffile:
        found = False
        infty = False
        for line in reffile:
            if '2 F(MAX):' in line:
                found = True
                words = line.split()
                if words[4] == 'Infinity':
                    infty = True

                    
    #checks if DCP was found
    if found and not infty:
        with open(f'trajectory-{trajectory}/DCP_{method}/fort.100') as reffile:
            no_DCP = False
            for line in reffile:
                if 'Phi:' in line:
                    words = line.split()
                    if words[1] == 'Infinity':
                        no_DCP = True
            
    #categorizes the calculated DCPs
    if found and not infty and not no_DCP:
        print(ev_deflection_check())
        R_new = reading_coordinates(trajectory, method)
        list_compared = compare_saddlepoints(R_new, trajectory) 

    elif not found:
        no_saddlepoint('no basin change')

    elif infty:
        no_saddlepoint('walked from basin to infty')

    elif no_DCP:
        no_saddlepoint('no DCP found')

#indices from list_compared: 0 = list_DCP, 1 = list_statistics, 2 = list_trajectories, 3 = list_enumerate_DCP, 4 = list_phi_values, 5 = list_order
enumeration = list_compared[3] + list_failure_DCP
statistics = list_compared[1] + list_failure_statistics
trajectories = list_compared[2] + list_failure_trajectories
Delta_phi = list_compared[4] + list_failure_phi_values
order = list_compared[5] + list_failure_order


#prints out .csv file
with open(f'DCP-analysis_{method}.csv', 'w') as reffile:
    reffile.write(f'runtime\t {runtime()}\t{runtime_hours(runtime())}\thours\t\n')
    reffile.write(f'\tfrequency\tDelta_phi\torder\ttrajectories\n')
    for i_item, item in enumerate(enumeration):
        reffile.write(f'{enumeration[i_item]}\t{statistics[i_item]}\t{Delta_phi[i_item]}\t{order[i_item]}\t{trajectories[i_item]}\n')

print(f'DCP-analysis_{method}.csv was generated')
