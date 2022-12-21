#!/bin/env python3

import numpy as np
from pyscript import *
import re
from sigfig import round

threshold = 1e-2


list_DCP = []
list_phi_values = []
list_enumerate_DCP =[]
list_statistics = []
list_trajectories = []
list_order = []

list_failure_DCP = []
list_failure_statistics = []
list_failure_trajectories = []
list_failure_phi_values = []
list_failure_order = []


def read_count():
    with open(f'trajectory.ami', 'r') as ami_file:
        notdone_count = True
        notdone_file = True
        for line in ami_file:
            if 'count' in line and notdone_count:
                count = int(re.search(r'\d+', line).group())
                notdone_count = False
    return(count)


def phi_value(x, trajectory):
    if x == 0:
        with open(f'./trajectory-start/minimum/fort.100') as reffile:
            found = False
            for line in reffile:
                found = True
                if 'Phi:' in line:
                   words = line.split()
                   phi = float(words[1])
            if not found:
                print('Phi from trajectory-start was not found')
    elif x == 1:
        with open(f'./trajectory-{trajectory}/DCP_{method}/fort.100') as reffile:
            found = False
            for line in reffile:
                if 'Phi:' in line:
                    words = line.split()
                    phi = float(words[1])
            if not found:
                print(f'Phi from trajectory-{trajectory}/DCP_{method} was not found')
    elif x == 2:
        with open(f'./trajectory-{trajectory}/minimum/fort.100') as reffile:
            found = False
            for line in reffile:
                if 'Phi:' in line:
                    found = True
                    words = line.split()
                    phi = float(words[1])
            if not found:
                print('Phi from trajectory-start was not found')
    return phi


def reading_order(trajectory):
    with open(f'trajectory-{trajectory}/DCP_{method}/fort.100') as reffile:
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


def compare_saddlepoints(R_new, trajectory):
    found = False
    for i_DCP, R_DCP in enumerate(list_DCP):
        norm = np.linalg.norm(R_new - R_DCP)
        if norm <= threshold:
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


#reads out and defines count and n_elecs
n_elecs = reading_n_elecs()
count = read_count()
            
#reads saddlepoint_calculation.in file (would be nicer here with regular expressions)
with open('saddlepoint_calculation.in', 'r') as reffile:
    for line in reffile:
        if 'threshold_DCP_guess' in line:
            words = line.split()
            threshold_DCP_guess = float(words[1])
        else:
            threshold_DCP_guess = 1e-1
            
        if 'method' in line:
            words = line.split()
            method = words[1]
        else: 
            method = 'gradient_norm'            
            
            
#starts evaluation
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
            
    if found and not infty and not no_DCP:
        with open(f'trajectory-{trajectory}/DCP_{method}/fort.100') as reffile:
            R_new = []
            for line in reffile:
                if 'after minimize:' in line:
                    line = reffile.readline()
                    line = reffile.readline()
                    for _ in range(n_elecs):
                        line = reffile.readline()
                        words = line.split()
                        for word in words[1:]:
                            R_new.append(float(word))

        compare_saddlepoints(np.array(R_new), trajectory)

    elif not found:
        found = False
        for i_failure, t_failure in enumerate(list_failure_DCP):
            if t_failure == 'no basin change':
                list_failure_statistics[i_failure] += 1
                list_failure_trajectories[i_failure] += f', {trajectory}'
                found = True
                break
        if not found:
            list_failure_DCP.append('no basin change')
            list_failure_statistics.append(1)
            list_failure_trajectories.append(f'{trajectory}')
            list_failure_phi_values.append('none')
            list_failure_order.append('none')

    elif infty:
        found = False
        for i_failure, t_failure in enumerate(list_failure_DCP):
            if t_failure == 'walked from basin to infty':
                list_failure_statistics[i_failure] += 1
                list_failure_trajectories[i_failure] += f', {trajectory}'
                found = True
                break
        if not found:
            list_failure_DCP.append('walked from basin to infty')
            list_failure_statistics.append(1)
            list_failure_trajectories.append(f'{trajectory}')
            list_failure_phi_values.append('none')
            list_failure_order.append('none')

    elif no_DCP:
        found = False
        for i_failure, t_failure in enumerate(list_failure_DCP):
            if t_failure == 'no DCP found':
                list_failure_statistics[i_failure] += 1
                list_failure_trajectories[i_failure] += f', {trajectory}'
                found = True
                break
        if not found:
            list_failure_DCP.append('no DCP found')
            list_failure_statistics.append(1)
            list_failure_trajectories.append(f'{trajectory}')
            list_failure_phi_values.append('none')
            list_failure_order.append('none')

enumeration = list_enumerate_DCP + list_failure_DCP
statistics = list_statistics + list_failure_statistics
trajectories = list_trajectories + list_failure_trajectories
Delta_phi = list_phi_values + list_failure_phi_values
order = list_order + list_failure_order



#prints out .csv file
with open(f'DCP-analysis_{method}.csv', 'w') as reffile:
    reffile.write(f'runtime\t {runtime()}\t{runtime_hours(runtime())}\thours\t\n')
    reffile.write(f'\tfrequency\tDelta_phi\torder\ttrajectories\n')
    for i_item, item in enumerate(enumeration):
        reffile.write(f'{enumeration[i_item]}\t{statistics[i_item]}\t{Delta_phi[i_item]}\t{order[i_item]}\t{trajectories[i_item]}\n')

print(f'DCP-analysis_{method}.csv was generated')
