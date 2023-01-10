#!/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

with open('saddlepoint_calculation.in', 'r') as reffile:
    found = False
    for line in reffile:  
        if 'method' in line:
            words = line.split()
            method = words[1]
            found = True
            break
    if not found:
        method = 'gradient_norm'       

#reads out the different DCP from DCP-analysis.csv
with open(f'DCP-analysis_{method}.csv', 'r') as reffile:
    list_barrier = []
    list_frequency = []
    list_label_temp = []
    list_frequency_temp = []
    list_adjacent = []
    line = reffile.readline()
    line = reffile.readline()
    for line in reffile:
        words = line.split('\t')
        list_barrier.append(words[2])
        list_frequency.append(int(words[1]))
        print(words[5])
        if words[5] == 'adjacent_minima':
            list_adjacent.append(True)
        else:
            list_adjacent.append(False)
        if 'none' in line:
            list_label_temp.append(words[0])
            list_frequency_temp.append(int(words[1]))            

print(list_adjacent)
            
#reads out the input information for the plot 
with open('assignment_graph.in', 'r') as reffile:
    list_category = []
    list_label = []
    for line in reffile:
        if 'category' in line:
            numbers = line.split()
            for number in numbers[1:]:
                list_category.append(number)
        if 'label' in line:
            words = line.split()
            for word in words[1:]:
                list_label.append(word)

#sums up the frequencies into the categories from the input file
list_category_frequency = []
list_category_frequency_adj = []
list_category_frequency_not_adj = []

for _ in list_category:
    list_category_frequency.append(0)
    list_category_frequency_adj.append(0)
    list_category_frequency_not_adj.append(0)
list_label.append('other_DCPs')
list_category_frequency.append(0)
list_category_frequency_adj.append(0)
list_category_frequency_not_adj.append(0)
    
for i_barrier, barrier in enumerate(list_barrier):
    if barrier not in list_category and barrier != 'none':
        list_category_frequency[-1] += list_frequency[i_barrier]
        if list_adjacent[i_barrier]:
            list_category_frequency_adj[i_category] += list_frequency[i_barrier]
        else:
            list_category_frequency_not_adj[i_category] += list_frequency[i_barrier]
    else:    
        for i_category, category in enumerate(list_category):
            if barrier == category:
                list_category_frequency[i_category] += list_frequency[i_barrier]
                if list_adjacent[i_barrier]:
                    list_category_frequency_adj[i_category] += list_frequency[i_barrier]
                else:
                    list_category_frequency_not_adj[i_category] += list_frequency[i_barrier]
                break
       
list_label += list_label_temp
list_category_frequency += list_frequency_temp
list_category_frequency_adj += list_frequency_temp
list_category_frequency_not_adj += list_frequency_temp
            
#sums up the barriers which are slightly different, but are set equal in the input file
list_final_label = []
list_final_frequencies = []
list_final_frequencies_adj = []
list_final_frequencies_not_adj = []
for i_label, label in enumerate(list_label):
    if label not in list_final_label:
        list_final_label.append(label)
        list_final_frequencies.append(list_category_frequency[i_label])
        list_final_frequencies_adj.append(list_category_frequency_adj[i_label])
        list_final_frequencies_not_adj.append(list_category_frequency_not_adj[i_label])
    else:
        list_final_frequencies[list_final_label.index(label)] += list_category_frequency[i_label]
        list_final_frequencies_adj[list_final_label.index(label)] += list_category_frequency_adj[i_label]
        list_final_frequencies_not_adj[list_final_label.index(label)] += list_category_frequency_adj[i_label]

#relabels import names from .csv
list_final_label[list_final_label.index('walked from basin to infty')] = 'to_infty'
list_final_label[list_final_label.index('no DCP found')] = 'no_DCP'
if 'no basin change' in list_final_label:
    list_final_label[list_final_label.index('no basin change')] = 'same_basin'
    
#creates output file
with open (f'assignment_graph_{method}.out', 'w') as printfile:
    printfile.write('categories: ')
    for item in list_final_label: 
        printfile.write(f'{item} ')
    printfile.write('\n')
    printfile.write('frequency: ')
    for item in list_final_frequencies: 
        printfile.write(f'{item} ')
    printfile.write('\n')
    printfile.write('frequency: adjacent')
    for item in list_final_frequencies_adj: 
        printfile.write(f'{item} ')
    printfile.write('\n')
    printfile.write('frequency: not adjacent')
    for item in list_final_frequencies_not_adj: 
        printfile.write(f'{item} ')
print(f'assignment_graph_{method}.out was generated') 

#makes bar diagram
#plt.rcParams['figure.figsize'] = (10, 8.5)
#plt.rcParams['font.size'] = 14
#plt.rcParams['lines.linewidth'] = 2
#plt.rc ('axes', titlesize = 20)
#plt.rc ('axes', labelsize = 18)
#plt.rc ('xtick', labelsize = 10) 
#
#plt.ylabel('frequency', labelpad = 15) 
#
#x_values = range(len(list_final_frequencies))
#
#plt.xticks(x_values, list_final_label)  
#
#plt.bar(x_values, list_final_frequencies, width = 0.5, label = list_label)
#
#plt.savefig(f'assignment_graph_{method}.png')
#print(f'assignment_graph_{method}.png was generated')
