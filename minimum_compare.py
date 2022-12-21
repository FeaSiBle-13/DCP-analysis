#!/bin/env python3
import matplotlib.pyplot as plt 
import numpy as np


def read_count():
    with open(f'trajectory.ami', 'r') as ami_file:
    notdone_count = True
    notdone_file = True
    for line in ami_file:
        if 'count' in line and notdone_count:
            count = int(re.search(r'\d+', line).group())
            notdone_count = False
    return(count)


def reading_n_elecs():
    with open(f'trajectory-1-max.ref', 'r') as reffile:
        for line in reffile:
            if '1 F(MAX):' in line:
                line = reffile.readline()
                words = line.split()
                n_elecs = int(words[0])
            else:
                print('n_elecs was not found in trajectory-1-max.ref')
    return(n_elecs)


def start_potential():
    with open(f'trajectory-1-max.ref', 'r') as newP_file:
        for line in reffile:
            if '1 F(MAX):' in line:
                words = newP_line.split()
    return(np.round(float(words[4]),5))


#reads in count (number of calculated trajectories), n_elecs and potential of start minimum
count = read_count()
n_elecs = reading_n_elecs()
start_potential = start_potential()

#reads out different potentials and electron movements for same potential as start potetnial and creates statistic
list_potentials = []
list_statistics = []
list_traj = []

list_elec_movements = []
list_statistics_mov = []
list_traj_mov = []

for trajectory in range(count+1):
    with open(f'trajectory-{trajectory}-max.ref','r') as reffile:
        line = reffile.readline()
        while f'MAX:' not in line:
            line = reffile.readline()
        line = reffile.readline()
        for _ in range(n_elecs+1):
            line = reffile.readline()    
    words = line.split()
    if len(words) < 5:
        print(f'WARNING: no second minimum in {trajectory} found')
    else:
        if np.round(float(words[4]),5) not in list_potentials:
            list_potentials.append(np.round(float(words[4]),5)) 
            list_statistics.append(1)
            list_traj.append(f'{trajectory}')
        else:
            list_statistics[list_potentials.index(np.round(float(words[4]),5))] +=1 
            list_traj[list_potentials.index(np.round(float(words[4]),5))] += f', {trajectory}'
        if round(float(words[4]),5) == round(start_potential, 5): 
            with open(f'trajectory-{trajectory}.idx') as idx_file:
                idx_line = idx_file.readline()
                idx_line = idx_file.readline()
                words_1 = idx_line.split()
                idx_line = idx_file.readline()
                words_2 = idx_line.split()
    
                elec_mov = ""
                for p in range(n_elecs):
                    if int(words_1[p+4]) != int(words_2[p+4]):
                        elec_mov += f'{int(words_1[p+4])} -> {int(words_2[p+4])}' + " " 
                if elec_mov not in list_elec_movements:
                    list_elec_movements.append(elec_mov)
                    list_statistics_mov.append(1)
                    list_traj_mov.append(f'{trajectory}')
                elif elec_mov in list_elec_movements:
                    list_statistics_mov[list_elec_movements.index(elec_mov)] +=  1
                    list_traj_mov[list_elec_movements.index(elec_mov)] += f', {trajectory}'


#prints out a .csv file
with open(f'minimum_compare.csv', 'w') as printfile:
    printfile.write(f'minimum potential\tDelta start potential\tmovement\tfrequency\ttrajectories\n')
    printfile.write(f'{start_potential}\t{int(0)}\t \t{int(0)}\tstart minimum\n')

    for i_item, item in enumerate(list_potentials):
        if item != start_potential:
            printfile.write(f'{item}\t{np.round(np.abs(start_potential-list_potentials[i_item]), 4)}\t \t{list_statistics[i_item]}\t{list_traj[i_item]}\n')
        else:
            printfile.write(f'{list_potentials[i_item]}\t{np.round(np.abs(start_potential-list_potentials[i_item]), 4)}\t \t{list_statistics[i_item]}\t{list_traj[i_item]}\n')
            for i_item, item in enumerate(list_elec_movements):
                printfile.write(f' \t \t{item}\t{int(list_statistics_mov[i_item])}\t{list_traj_mov[i_item]}\n')

print('minimum_compare.csv was generated')
