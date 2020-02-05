# Comparative-study-of-truss-sizing-optimisation-problems-in-MATLAB
Matlab codes of Comparative study of truss sizing optimisation problem 2020

Matlab codes of truss sizing optimisation including test problems, objective and constrained functions,
and demo of optimisation algorithms, SHAMODE and SHAMODE-WO are included.

Test problems are bi-objective truss sizing optimisation problems. Objective functions are mass and
compliance of structures calculated from finite element method. Design variables are cross-section areas
of each members of structures.

The codes can be used as benchmark test problems for developing of truss optimisation algorithms.

The optimisation run can be evaluated with 'main.m'. Output file of an optimisation run is in '.mat-file'
format e.g. 'rst_xxx_xxx.mat'. There are 4 metrics used to indicated performance of optimisation algorithms
included Hypervolume (HV), Generational Distance (GD), Inverted Generational Distance (IGD), and
Spacing-to-Extent (STE). Demo of metric calculations are provided in 'metric_calculation.m'.
