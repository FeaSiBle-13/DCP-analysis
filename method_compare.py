#!/bin/env python3

from pyscript import *
import re

with open(f'DCP-analysis_newton.csv', 'r') as reffile:
    line = reffile.readline()
    line = reffile.readline()
    for line in reffile:
        words = line.split('\t')
        trajectories = words[4].split(',')
        for item in trajectories:
            trajectory = re.search(r'\d+', item).group()
            print(f'traj = {trajectory}  pot = {words[1]}')
        
        
    
        
