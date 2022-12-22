#!/bin/env python3

from pyscript import *
import re

list_trajectory = []
list_potential = []

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
                        print(f'trajectory {item} changed from potential {list_potential[i_item]} to {potential}')
                        
                        
        
        
    
        
