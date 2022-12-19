#!/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

#reads out the different DCP from DCP-analysis.csv
with open(f'DCP-analysis.csv', 'r') as reffile:
    list_barrier = []
    list_frequency = []
    list_label_temp = []
    list_frequency_temp = []
    line = reffile.readline()
    line = reffile.readline()
    for line in reffile:
        words = line.split('\t')
        list_barrier.append(words[2])
        list_frequency.append(int(words[1]))
        if 'none' in line:
            list_label_temp.append(words[0])
            list_frequency_temp.append(int(words[1]))            

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

for _ in list_category:
    list_category_frequency.append(0)
list_label.append('other DCPs')
list_category_frequency.append(0)
    
for i_barrier, barrier in enumerate(list_barrier):
    for i_category, category in enumerate(list_category):
        if barrier == category:
            list_category_frequency[i_category] += list_frequency[i_barrier]
            break
        elif barrier not in list_category and barrier != 'none':
            list_category_frequency[-1] += list_frequency[i_barrier]

list_label += list_label_temp
list_category_frequency += list_frequency_temp            
            
#sums up the barriers which are slightly different, but are set equal in the input file
list_final_label = []
list_final_frequencies = []
for i_label, label in enumerate(list_label):
    if label not in list_final_label:
        list_final_label.append(label)
        list_final_frequencies.append(list_category_frequency[i_label])
    else:
        list_final_frequencies[list_final_label.index(label)] += list_category_frequency[i_label]
        
        


print(list_final_label)
print(list_final_frequencies) 
        

plt.rcParams['figure.figsize'] = (10, 8.5)
plt.rcParams['font.size'] = 14
plt.rcParams['lines.linewidth'] = 2
plt.rc ('axes', titlesize = 20)
plt.rc ('axes', labelsize = 18)
plt.rc ('xtick', labelsize = 10) 

#plt.title('Potentials $\Delta\phi$', pad= 20) 

#plt.xlabel('Maximum probability path', labelpad = 15) 
plt.ylabel('frequency', labelpad = 15) 

x_values = range(len(list_final_frequencies))

plt.xticks(x_values, list_final_label)  

#plt.xticks(np.arange(8, step=1), list_ordered_pot)  

plt.bar(x_values, list_final_frequencies, width = 1, label = list_label)

plt.savefig('saddlepoint_assignment_graph.png')
