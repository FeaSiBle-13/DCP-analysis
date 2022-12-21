#!/bin/env python3
import matplotlib.pyplot as plt 
import numpy as np

count = int(input('How many trajectories are to check (from 1- x) ?')) 

#reads out start potential and defines n_elecs
with open(f'trajectory-1-max.ref', 'r') as newP_file:
    newP_line = newP_file.readline()
    while f'MAX:' not in newP_line:
        newP_line = newP_file.readline()
    words = newP_line.split()
    p_start = np.round(float(words[4]),5)
    newP_line = newP_file.readline()
    words = newP_line.split()
    n_elecs = int(words[0])



#reads out different potentials and electron movements for same potential as start potetnial and creates statistic
list_potentials = []
list_statistics = []
list_traj = []

list_elec_movements = []
list_statistics_mov = []
list_traj_mov = []

for traj in range(count):
    with open(f'trajectory-{traj+1}-max.ref','r') as newP_file:
        newP_line = newP_file.readline()
        while f'MAX:' not in newP_line:
            newP_line = newP_file.readline()
        newP_line = newP_file.readline()
        for __ in range(n_elecs+1):
            newP_line = newP_file.readline()    
    newP_words = newP_line.split()
    if len(newP_words) < 5:
        print(f'WARNING: no second minimum in {traj + 1} found')
    else:
        if np.round(float(newP_words[4]),5) not in list_potentials:
            list_potentials.append(np.round(float(newP_words[4]),5)) 
            list_statistics.append(1)
            list_traj.append(f'{traj+1}')
        else:
            list_statistics[list_potentials.index(np.round(float(newP_words[4]),5))] +=1 
            list_traj[list_potentials.index(np.round(float(newP_words[4]),5))] += f', {traj+1}'
        if round(float(newP_words[4]),5) == round(p_start, 5): 
            with open(f'trajectory-{traj+1}.idx') as idx_file:
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
                    list_traj_mov.append(f'{traj+1}')
                elif elec_mov in list_elec_movements:
                    list_statistics_mov[list_elec_movements.index(elec_mov)] +=  1
                    list_traj_mov[list_elec_movements.index(elec_mov)] += f', {traj+1}'


#prints out a .csv file
with open(f'compare_combined.csv', 'w') as reffile:
    reffile.write(f'minimum potential\tDelta start potential\tmovement\tfrequency\ttrajectories\n')
    reffile.write(f'{p_start}\t{int(0)}\t \t{int(0)}\tstart minimum\n')

    for item in list_potentials:
        if item != p_start:
            reffile.write(f'{list_potentials[list_potentials.index(item)]}\t{np.round(np.abs(p_start-list_potentials[list_potentials.index(item)]), 4)}\t \t{list_statistics[list_potentials.index(item)]}\t{list_traj[list_potentials.index(item)]}\n')
        else:
            reffile.write(f'{list_potentials[list_potentials.index(item)]}\t{np.round(np.abs(p_start-list_potentials[list_potentials.index(item)]), 4)}\t \t{list_statistics[list_potentials.index(item)]}\t{list_traj[list_potentials.index(item)]}\n')
            for item_mov in list_elec_movements:
                reffile.write(f' \t \t{list_elec_movements[list_elec_movements.index(item_mov)]}\t{int(list_statistics_mov[list_elec_movements.index(item_mov)])}\t{list_traj_mov[list_elec_movements.index(item_mov)]}\n')


#makes a bad bar diagram
list_ordered_pot = []
list_ordered_stat = []

len_ = len(list_potentials)

for _ in range(len_):
    list_ordered_pot.append(min(list_potentials))
    list_ordered_stat.append(list_statistics[list_potentials.index(min(list_potentials))])
    if len(list_potentials) != 1:
        list_statistics.remove(list_statistics[list_potentials.index(min(list_potentials))])
        list_potentials.remove(min(list_potentials))
    

plt.rcParams['figure.figsize'] = (10, 8.5)
plt.rcParams['font.size'] = 14
plt.rcParams['lines.linewidth'] = 2 
plt.rc ('axes', titlesize = 20) 
plt.rc ('axes', labelsize = 18) 
plt.rc ('xtick', labelsize = 10) 


plt.title('Potentials $\Delta\phi$', pad= 20) 

plt.xlabel('Maximum probability path', labelpad = 15) 
plt.ylabel('$\Delta\phi$', labelpad = 15) 

iter_step = 0 
for freq in list_ordered_stat:
    plt.annotate(f'{freq}', (iter_step-0.1, freq+3))
    iter_step += 1

x_values = range(len(list_ordered_pot))
plt.xticks(np.arange(len(list_ordered_pot), step=1), list_ordered_pot)  

plt.bar(x_values, list_ordered_stat)
plt.savefig('compare_potentials.png')
