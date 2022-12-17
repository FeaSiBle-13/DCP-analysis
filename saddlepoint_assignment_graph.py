import numpy as np
import matplotlib.pyplot as plt

list_20 = [228, 6, 226, 171, 208, 2, 237, 1, 235, 211, 1, 7, 28, 1, 2, 1, 1, 1, 7, 26, 3, 3, 1, 4, 1, 1, 1, 4, 1, 1, 3, 2, 1, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 2, 4, 5, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 261, 13, 13]

list_barrier = [0.0337, 0.0326, 0.0337, 2.49e-06, 0.0337, 0.164, 0.0334, 2.04, 0.0337, 0.0334, 0.418, 0.0326, 0.164, 0.177, 0.0326, 9.32, 2.5, 0.177, 0.221, 0.164, 1.23, 0.0323, 1.51, 0.0328, 0.0438, 0.0438, 1.72, 0.0326, 12.1, 0.778, 0.0684, 0.0328, 1.16, 2.49e-06, 0.0893, 0.0678, 1.23, 2.04, 4.4, 1.23, 0.068, 0.236, 19.1, 1.23, 0.0328, 0.0323, 5.88, 0.0678, 0.434, 1.43, 3.44, 0.182, 0.0328, 3.32, 1.23, 2.5, 0.0326, 1.23, 0.23, 0.068, 2.2, 0.0667, 2.04, 0.0328, 2.04, 5.58, 1.23, 25.6, 0.0684, 0.0326, 0.0323, 10.9, 11.5, 17.9, 1.28, 0.181, 1.99, 3.16, 0.068, 0.068, 2.04, 4.04, 1.84, 0.759, 8.84, 3.92, 0.0932, 1.72, 1.51, 1.18, 2.04, 1.52, 0.0138, 1.87, 0.386, 1.36, 4.18, 28.9, 2.04, 1.51, 'none', 'none', 'none']


CH_ion = 0
three_same_spin = 0
CC_cov = 0
not_dedicated = 0


for i_barrier, barrier in enumerate(list_barrier):
    if barrier == 0.0337 or barrier == 0.0334:
        CH_ion += list_20[i_barrier]
    elif barrier == 1.23:
        three_same_spin += list_20[i_barrier]
    elif barrier == 0.164:
        CC_cov += list_20[i_barrier]
    elif barrier != 0.0337 or barrier != 0.0334 or barrier != 1.23 or barrier != 0.164 or barrier != 'none':
        not_dedicated += list_20[i_barrier] 
    
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
