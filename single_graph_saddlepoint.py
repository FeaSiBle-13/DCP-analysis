#!/bin/env python3

count = int(input("Which trajectory do you mean?"))
import matplotlib.pyplot as plt #module load matplotlib.pyplot
from sigfig import round #pip install sigfig
import numpy as np

with open(f'./trajectory-start/minimum/fort.100') as reffile:
    line= reffile.readline()
    while 'Phi:' not in line:
        line = reffile.readline()
    words = line.split()
    phi_SCP1 = float(words[1])


with open(f'./trajectory-{count}/DCP_newton/fort.100') as pot1_infile:
    pot1_line =  pot1_infile.readline()
    while 'Phi:' not in pot1_line:
        pot1_line = pot1_infile.readline()
    pot1_words= pot1_line.split()
    phi_DCP = float(pot1_words[1])
    Dphi_DCP= np.abs(phi_SCP1-phi_DCP)


with open(f'./trajectory-{count}/minimum/fort.100') as pot2_infile:
    pot2_line =  pot2_infile.readline()
    while 'Phi:' not in pot2_line:
        pot2_line = pot2_infile.readline()
    pot2_words= pot2_line.split()
    phi_SCP2 = float(pot2_words[1])
    Dphi_SCP2= np.abs(phi_SCP1-phi_SCP2)

norm1= min(phi_SCP1, phi_DCP, phi_SCP2)
norm2= max(phi_SCP1, phi_DCP, phi_SCP2)
phi_SCP1_norm = (phi_SCP1-norm1)/(norm2-norm1)
phi_DCP_norm = (phi_DCP-norm1)/(norm2-norm1)
phi_SCP2_norm = (phi_SCP2-norm1)/(norm2-norm1)

x_values = [-0.05, 0.05, 0.95, 1.05, 1.95, 2.05]
phi_values_normed = [phi_SCP1_norm, phi_SCP1_norm, phi_DCP_norm, phi_DCP_norm, phi_SCP2_norm, phi_SCP2_norm]

plt.rcParams['figure.figsize'] = (10, 8.5)
plt.rcParams['font.size'] = 14
plt.rcParams['lines.linewidth'] = 2
plt.rc('axes', titlesize = 20)
plt.rc('axes', labelsize = 18)

plt.title('Potentials $\Delta\phi$', pad= 20)

plt.xlabel('Maximum probability path', labelpad = 15)
plt.ylabel('$\Delta\phi$', labelpad = 15)

plt.annotate(f'$\Delta\phi$ = {round(Dphi_DCP, sigfigs=3)}',(0.8, phi_DCP_norm +0.12))
plt.annotate(f'$\phi$ = {round(phi_DCP, sigfigs = 6)}',(0.84, phi_DCP_norm + 0.06))
plt.annotate(f'$\Delta\phi$ = {round(Dphi_SCP2, sigfigs = 3)}',(1.72, phi_SCP2_norm-0.06))
plt.annotate(f'$\phi$ = {round(phi_SCP2, sigfigs = 6)}',(1.76, phi_SCP2_norm-0.12))
plt.annotate(f'$\phi$ = {round(phi_SCP1, sigfigs = 6)}',(-0.08, phi_SCP1_norm-0.06))

ax = plt.gca()
plt.ylim(-0.19, 1.2)
#plt.yticks()
ax.axes.yaxis.set_ticks([0])

labeling = ['starting SCP', 'DCP', 'ending SCP']
plt.xticks(np.arange(0, 2.25, step=1), labeling)


plt.plot(x_values[1:3], phi_values_normed[1:3], color='blue', linestyle='--', marker='')
plt.plot(x_values[3:5], phi_values_normed[3:5], color='blue', linestyle='--', marker='')
plt.plot(x_values[2:4], phi_values_normed[2:4], color = 'blue')
plt.plot(x_values[0:2], phi_values_normed[0:2], color = 'blue')
plt.plot(x_values[4:7], phi_values_normed[4:7], color = 'blue')

#plt.show()
plt.savefig(f'trajectory-{count}/result/potentialbarrier_traj-{count}.png')
