#!/bin/env python3

from pyscript import *

with open(f'DCP-analysis_newton.csv', 'r') as reffile:
    line = reffile.readline()
    line = reffile.readline()
    for line in reffile:
        words = line.split('\t')
        trajectories = words[4].split(',')
        for trajectory in trajectories:
            print(f'traj = {trajectory}  pot = {words[1]}')
        
        
    
        
