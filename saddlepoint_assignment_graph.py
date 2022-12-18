import numpy as np
import matplotlib.pyplot as plt

#reads out the different DCP from DCP-analysis.csv
with open(f'DCP-analysis.csv', 'r') as reffile:
    list_barrier = []
    list_frequency = []
    line = reffile.readline()
    line = reffile.readline()
    for line in reffile:
        words = line.split()
        list_barrier.append(words[2])
        list_frequency.append(words[1])
        

CH_ion = 0
three_same_spin = 0
CC_cov = 0
not_dedicated = 0


for i_barrier, barrier in enumerate(list_barrier):
    if barrier == 0.0337 or barrier == 0.0334:
        CH_ion += list_frequency[i_barrier]
    elif barrier == 1.23:
        three_same_spin += list_frequency[i_barrier]
    elif barrier == 0.164:
        CC_cov += list_frequency[i_barrier]
    elif barrier != 0.0337 or barrier != 0.0334 or barrier != 1.23 or barrier != 0.164 or barrier != 'none':
        not_dedicated += list_frequecy[i_barrier] 
    
DCP_notfound = 261
no_basin_change = 13
walked_to_infty = 13

        
list_result = [CH_ion, three_same_spin, CC_cov, not_dedicated-(DCP_notfound + walked_to_infty + no_basin_change), DCP_notfound, no_basin_change, walked_to_infty ] 
list_label = ['CH_ion', 'three_same_spin', 'CC_cov', 'not_dedicated', 'DCP not found', 'no basin change', 'walked from basin to infty'] 
    
print(list_result)
print(list_label)
    


print(f'sum = {CH_ion + three_same_spin + CC_cov + not_dedicated}')
        
list_20 = [1312, 30, 59, 269, 316, 0, 14]
list_17 = [1342, 24, 68, 262, 283, 0, 21]
list_15 = [1345, 14, 56, 298, 261, 13, 13]
list_label = ['CH_ion', '3 same', 'CC_cov', 'not_ddctd', 'no DCP', 'same basin', 'to infty']
    


plt.rcParams['figure.figsize'] = (10, 8.5)
plt.rcParams['font.size'] = 14
plt.rcParams['lines.linewidth'] = 2
plt.rc ('axes', titlesize = 20)
plt.rc ('axes', labelsize = 18)
plt.rc ('xtick', labelsize = 10) 

#plt.title('Potentials $\Delta\phi$', pad= 20) 

 

#plt.xlabel('Maximum probability path', labelpad = 15) 
plt.ylabel('frequency', labelpad = 15) 

x_values = range(len(list_20))

x_values15 = []
for x in range(len(list_15)):
    x_values15.append(x-0.25)

x_values20 = []
for x in range(len(list_20)):
    x_values20.append(x+0.25)

plt.xticks(x_values, list_label)  



#plt.xticks(np.arange(8, step=1), list_ordered_pot)  


plt.bar(x_values20, list_20, width = 0.25, label = legend_con)
plt.bar(x_values, list_17, width = 0.25)
plt.bar(x_values15, list_15, width = 0.25)

plt.savefig('bar_Ethane_stat-22-12-15.png')
