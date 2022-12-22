#!/bin/env python3

from pyscript import *

with open(f'DCP-analysis_newton.csv', 'r') as reffile:
    line = reffile.readline()
    line = reffile.readline()
    for line in reffile:
        words = line.split('\t')
        trajectories = int(words[4].split(,))
        print(f'traj = trajectories  pot = {words[1]}')
        
    
        
