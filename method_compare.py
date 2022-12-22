from pyscript import *

with open(f'DCP-analysis_newton', 'r') as reffile:
    line = reffile.readline()
    line = reffile.readline()
    for line in reffile:
        words = line.split()
        print(f'traj = {words[3]}  pot = {words[1]}')
        
    
        
