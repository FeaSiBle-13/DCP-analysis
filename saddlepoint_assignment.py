#!/bin/env python3
import subprocess

import numpy as np
from pyscript import *
import re
from sigfig import round

list_failure_DCP = []
list_failure_statistics = []
list_failure_trajectories = []
list_failure_phi_values = []
list_failure_order = []
list_failure_DCP_location = []


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


def phi_value(point_at_MPP, trajectory):
    if point_at_MPP == 'start':
        path = f'trajectory-start/minimum/fort.100'
    elif point_at_MPP == 'DCP':
        path = f'trajectory-{trajectory}/DCP_{method}/fort.100'
    elif point_at_MPP == 'end':
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


def read_order(trajectory):
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


def read_n_elecs():
    with open(f'trajectory-1-max.ref', 'r') as reffile:
        for line in reffile:
            if 'MAX:' in line:  
                line = reffile.readline()
                words = line.split()
                n_elecs = int(words[0])
    return n_elecs


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
    
    
def read_eigenvector(trajectory):
    if method == 'gradient_norm':
        path = f'trajectory-{trajectory}/DCP_{method}/newton_singlepoint/fort.100'
    if method == 'newton':
        path = f'trajectory-{trajectory}/DCP_newton/fort.100'
    with open(path) as reffile:
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
    cp(f'{name}.wf', f'eigenvector_check')
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


def compare_position(R1, R2, compare_mode):
    norm = np.linalg.norm(R1 - R2)
    if compare_mode == 'molecule_wise':
        if norm <= threshold_molecule:
            return(True)
        else:
            return(False)
    elif compare_mode == 'electron_wise':
        for l in range(n_elecs):
            same = True
            norm = np.linalg.norm(R1[l*3:l*3+3]-R2[l*3:l*3+3])
            if norm > threshold_electrons:
                same = False
                break
                return(False)
        if same:
            return(True)
    
    
def ev_deflection_check():
    eigenvec = read_eigenvector(trajectory)
    saddlepoint = read_coordinates(trajectory, method)
    
    #compares if stepest descent minimization ends in the starting or ending minimum
    list_minimized_deflection = []
    
    for m in range(1, 3):
        stepest_descent(trajectory, name, deflection_saddlepoint(eigenvec, saddlepoint, (-1) ** m * deflection_factor))
        list_minimized_deflection.append(read_coordinates(trajectory, 'stedes_eigvec'))
        for obj in ls(f'eigenvector_check'):
            rm(f'eigenvector_check/{obj}')
    
    list_found_min = [False, False]
    for m in range(1, 3):
        for vectors in list_minimized_deflection:
            same = compare_position(read_ref_file('max', m), vectors,  compare_mode)
            if same:
                list_found_min[m-1] = True
                break
    
    if list_found_min[0] and list_found_min[1]:
        return('adjacent_minima')
    else:
        return('no_adjacent_minima')
    

def compare_saddlepoints(R_new, trajectory, list_compared, ev_deflection_check):
    #indices from list_compared: 0 = list_DCP, 1 = list_statistics, 2 = list_trajectories, 3 = list_enumerate_DCP, 4 = list_phi_values, 5 = list_order, 6 = saddlepoint_location 
    found = False
    for i_DCP, R_DCP in enumerate(list_compared[0]):
        same = compare_position(R_new, R_DCP, compare_mode)
        if same and ev_deflection_check == list_compared[6][i_DCP]:
            list_compared[1][i_DCP] += 1
            list_compared[2][i_DCP] += f', {trajectory}'
            found = True
            break
    if not found:
        list_compared[0].append(R_new)
        list_compared[1].append(1)
        list_compared[2].append(f'{trajectory}')
        list_compared[3].append(len(list_compared[3])+1)
        list_compared[4].append(round(abs(phi_value('start', trajectory)-phi_value('DCP', trajectory)), sigfigs = 3 ))
        list_compared[5].append(f'{read_order(trajectory)}. order')
        list_compared[6].append(ev_deflection_check)
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
        list_failure_DCP_location.append('none')
        

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


#script starts here
#reads out and defines count and n_elecs
n_elecs = read_n_elecs()
count = read_trajectory_ami('count')
name = read_trajectory_ami('file')
            
#reads DCP_analysis.in file (would be nicer here with regular expressions)
threshold_molecule = read_DCP_analysis_in('threshold_molecule')
method = read_DCP_analysis_in('method')
deflection_factor = read_DCP_analysis_in('deflection_factor')
threshold_electrons = read_DCP_analysis_in('threshold_electrons')
compare_mode = read_DCP_analysis_in('compare_mode')


            
#starts evaluation
try:
    rm('eigenvector_check', True)
except subprocess.CalledProcessError:
    pass

mkdir('eigenvector_check')
list_compared = [[], [], [], [], [], [], []]

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
                if 'converged=' in line:
                    words = line.split()
                    if words[1] == 'F':
                        no_DCP = True
            
    #categorizes the calculated DCPs
    list_compared_old = list_compared
    if found and not infty and not no_DCP:
        R_new = read_coordinates(trajectory, method)
        list_compared = compare_saddlepoints(R_new, trajectory, list_compared_old, ev_deflection_check()) 
    elif not found:
        no_saddlepoint('no basin change')

    elif infty:
        no_saddlepoint('walked from basin to infty')

    elif no_DCP:
        no_saddlepoint('no DCP found')

#indices from list_compared: 0 = list_DCP, 1 = list_statistics, 2 = list_trajectories, 3 = list_enumerate_DCP, 4 = list_phi_values, 5 = list_order 6 = saddlepoint_location
enumeration = list_compared[3] + list_failure_DCP
statistics = list_compared[1] + list_failure_statistics
trajectories = list_compared[2] + list_failure_trajectories
Delta_phi = list_compared[4] + list_failure_phi_values
order = list_compared[5] + list_failure_order
DCP_location = list_compared[6] + list_failure_DCP_location

rm('eigenvector_check', True)

#frequency of saddlepoints which are between the minima from *max.ref file
frequency_adjacent = 0
for i_item, item in enumerate(DCP_location):
    if item == 'adjacent_minima':
        frequency_adjacent += statistics[i_item]
        

#prints out .csv file
with open(f'saddlepoint_assignment_{method}.csv', 'w') as reffile:
    reffile.write(f'runtime\t {runtime()}\t{runtime_hours(runtime())}\thours\t\t\n')
    reffile.write(f'\tfrequency\tDelta_phi\torder\ttrajectories\tsaddlepoint_location\n')
    for i_item, item in enumerate(enumeration):
        reffile.write(f'{enumeration[i_item]}\t{statistics[i_item]}\t{Delta_phi[i_item]}\t{order[i_item]}\t{trajectories[i_item]}\t{DCP_location[i_item]}\n')
    reffile.write('\n')
    reffile.write(f'frequency_adjacent\t{frequency_adjacent}\n')
print(f'saddlepoint_assignment_{method}.csv was generated')
