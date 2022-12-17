#!/bin/env python3

import subprocess

from pyscript import *
import re
import numpy as np

threshold_DCP_guess = 1e-1


def reading_n_elecs():
    with open(f'trajectory-1-max.ref', 'r') as reffile:
        line = reffile.readline()
        while 'MAX:' not in line:
            line = reffile.readline()
        line = reffile.readline()
        words = line.split()
        n_elecs = int(words[0])
    return n_elecs


def starting_maximum(trajectory):
    with open(f'trajectory-{trajectory}-max.ref', 'r') as reffile:
        R_max = []
        line = reffile.readline()
        while '1 F(MAX):' not in line:
            line = reffile.readline()
        line = reffile.readline()
        for _ in range(n_elecs):
            line = reffile.readline()
            words = line.split()
            for word in words:
                R_max.append(float(word))
        return np.array(R_max)


def ending_maximum(trajectory):
    with open(f'trajectory-{trajectory}-max.ref', 'r') as reffile:
        R_max = []
        line = reffile.readline()
while '2 F(MAX):' not in line:
            line = reffile.readline()
        line = reffile.readline()
        for _ in range(n_elecs):
            line = reffile.readline()
            words = line.split()
            for word in words:
                R_max.append(float(word))
        return np.array(R_max)
        print(R_max)


def basin_left(trajectory):
    with open(f'./trajectory-{trajectory}-traj.ref', 'r') as reffile:
        R_i = []
        line = reffile.readline()
        while '1 F(TRAJ):' not in line:
            line = reffile.readline()
        line = reffile.readline()
        for _ in range(n_elecs):
            line = reffile.readline()
            words = line.split()
            for word in words:
                R_i.append(float(word))
    return np.array(R_i)


def basin_enter(trajectory):
    with open(f'./trajectory-{trajectory}-traj.ref', 'r') as reffile:
        R_j = []
        line = reffile.readline()
        while '2 F(TRAJ):' not in line:
            line = reffile.readline()
        line = reffile.readline()
        for _ in range(n_elecs):
            line = reffile.readline()
            words = line.split()
            for word in words:
                R_j.append(float(word))
    return np.array(R_j)
 
 def DCP_guess_point(R_max1, R_max2, threshold_DCP_guess):
    for l in range(n_elecs):
        norm = np.linalg.norm(R_max1[l*3:l*3+3]-R_max2[l*3:l*3+3])
        if norm < threshold_DCP_guess:
            for r in R_max2[l*3:l*3+3]:
                printfile.write(f'{r} ')
        else:
            for t in R_int[l*3:l*3+3]:
                printfile.write(f'{t} ')
        printfile.write('\n')


#reads count from .ami
with open(f'trajectory.ami', 'r') as ami_file:
    notdone_count = True
    notdone_file = True
    for line in ami_file:
        if 'count' in line and notdone_count:
            count = int(re.search(r'\d+', line).group())
            notdone_count = False
        elif 'file' in line and notdone_file:
            name = re.search(r'file=([\'"]?)(.+?)\.wf\1', line).group(2)
            notdone_file = False

for trajectory in range(1, count+1):
    with open(f'trajectory-{trajectory}-max.ref') as reffile:
        found = False
        infty = False
        for line in reffile:
            if '2 F(MAX):' in line:
                found = True
                words = line.split()
                if words[4] == 'Infinity':
                    infty = True

    if found and not infty:
        #calculates the minima
        #makes folders
        mkdir(f'trajectory-{trajectory}')
        mkdir(f'trajectory-{trajectory}/result')
        cp(f'trajectory-start/result/cluster_start-out.yml', f'trajectory-{trajectory}/result')
        mkdir(f'trajectory-{trajectory}/minimum')
        for obj in ls('trajectory-source/minimum'):
            cp(f'trajectory-source/minimum/{obj}', f'trajectory-{trajectory}/minimum')

        #reads out the coordinates of the starting minimum
        n_elecs = reading_n_elecs()
        R_max = ending_maximum(trajectory)

        #creates ami file for Amolqc
        with open(f'trajectory-{trajectory}/minimum/stedes.ami', 'w') as sted_file:
            sted_file.write(f'''! seed for random number generation, not important
$gen(seed=101)
! reading the wave function file
$wf(read,file='{name}.wf')
$init_rawdata_generation()
$init_max_search(
step_size=0.1,
correction_mode=cut,
singularity_threshold=0.001,
method=steepest_descent,
verbose=2,
negative_eigenvalues=0,
eigenvalue_threshold=1e-10)
! setting the initial position
$init_walker(
free
''')
            for s in range(n_elecs):
                for t in R_max[s*3:s*3+3]:
                    sted_file.write(f'{t} ')
                sted_file.write('\n')
            sted_file.write(''')
$sample(create, size=1, single_point)
! maximize the walker
$maximize_sample()''')

        #makes the run
        with cd(f'trajectory-{trajectory}/minimum'):
            run('amolqc stedes.ami')
            run('ProcessMaxima cluster.yml')
            cp('cluster-out.yml', '../result/cluster_end-out.yml')



#calculates saddlepoints with newton
        #makes folders
        mkdir(f'trajectory-{trajectory}/DCP_newton')
        for obj in ls('trajectory-source/DCP_newton'):
            cp(f'trajectory-source/DCP_newton/{obj}', f'trajectory-{trajectory}/DCP_newton')


        #interpolation of both -traj- vectors
        R_int = (basin_left(trajectory)+basin_enter(trajectory))/2


        #creates the newton.ami with the interpolated vectors and the input with the fixed e positions
        with open(f'trajectory-{trajectory}/DCP_newton/newton.ami', 'w') as printfile:
            printfile.write(f'''! seed for random number generation, not important
$gen(seed=101)
! reading the wave function file
$wf(read,file='{name}.wf')
$init_rawdata_generation()
$init_max_search(
step_size=0.1,
correction_mode=none,
singularity_threshold=0.0001,
method=newton,
verbose=2,
negative_eigenvalues=-1,
eigenvalue_threshold=1e-10)
! setting the initial position
$init_walker(
free
''')
            DCP_guess_point(starting_maximum(trajectory), ending_maximum(trajectory), threshold_DCP_guess)
            printfile.write(''')
$sample(create, size=1, single_point)
! maximize the walker
$maximize_sample()''')

#executes amolqc and Processmaxima
        with cd(f'trajectory-{trajectory}/DCP_newton'):
            success = True
            try:
                run('amolqc newton.ami')
            except subprocess.CalledProcessError:
                success = False

            if success:
                run('ProcessMaxima cluster.yml')
                cp('cluster-out.yml', '../result/cluster_trans-out.yml')
            else:
                print(f'newton method did not work for trajectory-{trajectory}')


        #executesEVanalysis
        if success:
            with cd(f'trajectory-{trajectory}/result'):
                try:
                    run(f'/home/theochem/Scripts/EVanalysis/EVanalysis.py ../minimum/{name}.wf cluster_trans-out.yml')
                except subprocess.CalledProcessError:
                    pass


    else:
        print(f'trajectory-{trajectory} could not be calculated')                                                     


