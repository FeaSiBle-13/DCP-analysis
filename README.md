# DCP-analysis

### saddlepoint_calculation
After a trajectory run with Amolqc and after a run with minimum_calcuation.py, saddlepoints can be calculated with this script. An input file 'saddlepoint_calculation.in' is required. For example:

> threshold_DCP_guess= 1e-1  
  method= gradient_norm 
  
'threshold_DCP_guess' defines above which value the electron positions of two vectors are set as different.
'method' can be 'newton', 'gradient_norm' or 'none'.

### saddlepoint_assignment_graph

An assignment_graph.in file is required. It needs to contain a line with 'category:' followed by the barriers which should be shown as bars in the bar diagram. As seperator a space is used. A line with 'label:' is required where the categories are named. If several categories are accorded to one bar (so the same label) just label the categories in line label with the same name.

