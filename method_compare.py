#!/bin/env python3

from pyscript import *
import re

list_trajectory = []
list_potential = []
list_trajectory_change = []
list_potential_gradient = []
list_potential_newton = []


def print_statistic(list_potential_gradient, list_potential_newton):
    list_statistic = []
    list_stat_potential_gradient = []
    list_stat_potential_netwon = []
    for i_item, item in enumerate(list_potential_newton):
        for i_element, element in enumerate(list_stat_potential_newton):
            if item == element and list_potential_newton[item] == list_ 
            


with open(f'DCP-analysis_newton.csv', 'r') as reffile:
    line = reffile.readline()
    line = reffile.readline()
    for line in reffile:
        words = line.split('\t')
        trajectories = words[4].split(',')
        for item in trajectories:
            trajectory = re.search(r'\d+', item).group()
            potential = words[2]
            list_trajectory.append(trajectory)
            list_potential.append(potential)

with open(f'DCP-analysis_gradient_norm.csv', 'r') as reffile:
    line = reffile.readline()
    line = reffile.readline()
    for line in reffile:
        words = line.split('\t')
        trajectories = words[4].split(',')
        for item in trajectories:
            trajectory = re.search(r'\d+', item).group()
            potential = words[2]
            for i_item, item in enumerate(list_trajectory):
                if item == trajectory:
                    if potential != list_potential[i_item]:
                        list_trajectory_change.append(item)
                        list_potential_newton.append(list_potential[i_item])
                        list_potential_gradient.append(potential)

 
    
with open('method_compare.out', 'w') as printfile:
    printfile.write('trajectory\t potential_newton\t potential_gradient_norm\n')
    for i_item, item in enumerate(list_trajectory_change):
        printfile.write(f'{item}\t{list_potential_newton[i_item]}\t{list_potential_gradient[i_item]}\n')
  
print('file method_compare.out was generated')
                       
                        
                        
        
        
    
        
