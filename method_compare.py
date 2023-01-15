#!/bin/env python3

from pyscript import *
import re

list_trajectory = []
list_potential = []
list_trajectory_change = []
list_potential_gradient = []
list_potential_newton = []
list_adjacent_newton = []
list_adjacent_gradient = []
list_stat_adjacent_newton = []
list_stat_adjacent_gradient = []

with open(f'DCP-analysis_newton.csv', 'r') as reffile:
    line = reffile.readline()
    line = reffile.readline()
    for line in reffile:
        words = line.split('\t')
        if len(words) >= 5: 
            trajectories = words[4].split(',')
            for item in trajectories:
                trajectory = re.search(r'\d+', item).group()
                potential = words[2]
                searched = words[5]
                list_trajectory.append(trajectory)
                list_potential.append(potential)
                list_adjacent_newton.append(searched)

with open(f'DCP-analysis_gradient_norm.csv', 'r') as reffile:
    line = reffile.readline()
    line = reffile.readline()
    for line in reffile:
        words = line.split('\t')
        if len(words) >= 5:
            trajectories = words[4].split(',')
            for item in trajectories:
                trajectory = re.search(r'\d+', item).group()
                potential = words[2]
                searched = words[5]
                for i_item, item in enumerate(list_trajectory):
                    if item == trajectory:
                        if potential != list_potential[i_item]:
                            list_trajectory_change.append(item)
                            list_potential_newton.append(list_potential[i_item])
                            list_potential_gradient.append(potential)
                            list_adjacent_newton.append(searched)
                        

list_statistic = []
list_stat_potential_gradient = []
list_stat_potential_newton = []
for i_item, item in enumerate(list_potential_newton):
    empty = False
    if len(list_statistic) == 0:
        list_stat_potential_newton.append(item)
        list_stat_potential_gradient.append(list_potential_gradient[i_item])
        list_statistic.append(1)
        list_stat_adjacent_newton.append(list_adjacent_newton[i_item])
        list_stat_adjacent_gradient.append(list_adjacent_gradient[i_item])
        empty = True
    if not empty:
        for i_element, element in enumerate(list_stat_potential_newton):
            found = False
            if item == element and list_potential_gradient[i_item] == list_stat_potential_gradient[i_element] and list_adjacent_newton[i_item] == list_stat_adjacent_newton[i_element] and list_adjacent_gradient[i_item] == list_stat_adjacent_gradient[i_element] :
                list_statistic[i_element] += 1
                found = True
                break
        if not found:
            list_stat_potential_newton.append(item)
            list_stat_potential_gradient.append(list_potential_gradient[i_item])
            list_statistic.append(1)
            list_stat_adjacent_newton.append(list_adjacent_newton[i_item])
            list_stat_adjacent_gradient.append(list_adjacent_gradient[i_item])
    
    
with open('method_compare.out', 'w') as printfile:
    printfile.write('trajectory\t potential_newton\t potential_gradient_norm\t adjacent newt\tadjacent grad\tstatistic\n')
    for i_item, item in enumerate(list_statistic):
        printfile.write(f'\t{list_stat_potential_newton[i_item]}\t{list_stat_potential_gradient[i_item]}\t{list_stat_adjacent_newton}\t {list_stat_adjacent_gradient}\t{item}\n')
    printfile.write('\t \t \t \t \t \n')
    for i_item, item in enumerate(list_trajectory_change):
        printfile.write(f'{item}\t{list_potential_newton[i_item]}\t{list_potential_gradient[i_item]}\t{list_stat_adjacent_newton}\t{list_stat_adjacent_gradient}\n')
  
print(list_statistic)
print(list_stat_potential_newton)
print(list_stat_potential_gradient)
print(list_stat_adjacent_newton)
print(list_stat_adjacent_gradient)
print('file method_compare.out was generated')
                       
                        
                        
        
        
    
        
