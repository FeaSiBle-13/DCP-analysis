# DCP-analysis

## minimum_calculation
After a trajectory run with Amolqc, the second minima of a trajectory run are calculated with this script. The script is to be started in the folder where the Amolqc run was started. The results are to find in the folder trajectory-X/result. 


### saddlepoint_calculation
After a trajectory run with Amolqc and after a run with minimum_calcuation.py, saddlepoints can be calculated with this script. The script is to be started in the folder where the Amolqc run was started. An input file 'saddlepoint_calculation.in' is required. For example:

> threshold_DCP_guess= 1e-1  
  method= gradient_norm 
  
'threshold_DCP_guess' defines above which value the electron positions of two vectors are set as different.
'method' can be 'newton', 'gradient_norm' or 'none'.
As seperator a space is used after '='.


### saddlepoint_assignment
After a run with saddlepoint_calculation.py the saddlepoints are compared which eachother and are summarized in groups of same saddlepoints with this script. The script is to be started in the folder where the Amolqc run was started. A DCP-analysis_method.csv file is generated, in which the saddlepoints are summarized with information about the **frequency**, the **potentialbarrier**, **order** and the **trajectories** according to the group. Additionally the **runtime** for the saddlepoint_calculation run is given.

### saddlepoint_categorization
The saddlepoint_categorization.py script can be used, when a DCP-analysis_method.csv file was generated with the script saddlepoint_assignment. The script is to be started in the folder where the Amolqc run was started.
An assignment_graph.in file is required. For example:

> category: 0.0337 0.0334 0.164 1.23 0.0326 0.0328 0.0323  
  label: CH_ion CH_ion CC_cov 3_same_spin CH_ion CH_ion CH_ion
 
As seperator a space is used. 
It needs to contain a line with 'category:' and a line with 'label:'. After 'category' the categories are defined in which the saddlepoints should be summarized in. With 'label' the category names are defined. If several categories are accorded to the same label  the categories need to have the same label name. 
A assignment_graph.out file is generated with **categories** and the **frequency** of the saddlepoints according to the categories. The category **other_DCPs** is always given, when potentials are found, which are not given in the input file as 'category'. The categories **no_DCP**, **to_infty** and **same_basin** are generated, when the potential is 'none' in the DCP-analysis_method.csv file. The .out file is for example:

>categories: CH_ion CC_cov 3_same_spin other_DCPs no_DCP to_infty  
 frequency: 1423 45 15 495 8 14

Additionally an assignment_graph_method.png file with a bar diagram is generated.
