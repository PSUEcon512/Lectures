
CODE README FILE:  MPEC_LOG (National Border Only)

This file provides information about the the code release which accompanies "BORDERS, GEOGRAPHY AND OLIGOPOLY: EVIDENCE FROM THE WIND TURBINE INDUSTRY" by Kerem Cosar, Paul Grieco, and Felix Tintelnot.

The primary code was run uses MATLAB version 2012b and KNITRO 8.1.1 and 8.2. It has been run successfully on a Mac Pro OSX version 8.10.5 and on the Lion-XF cluster at Pennsylvania State University (https://rcc.its.psu.edu/resources/hpc/lionxf/#specs).  The bootstrap computation is parallelized and makes use of the PBS batch system on the Lion-XF cluster. 

What follows is a description of the files in this directory, which estimates the specification with only a national border: 

CODE FILES: 
Code files are listed in order in which they are called: 

master.m - Master file which runs code in order of dependencies. Running this file will run the specification and save the result.

setup.m - This file reads in the the data files and constructs the structures which will be used in estimation.  All model structures are saved in the output file "data_structures.mat" 

Main_mpec.m - Constructs sparsity pattern, calls optimization, and computes standard errors for the primary estimation. Saves results in 'est_point.mat'

dummy_objective.m - A dummy function which returns zero.  Used when we want to simply
solve for equilibrium instead of optimize the likelihood.  It is used as an objective
function input to KNITRO where the equilibrium constraints are supplied as constraints. 

model_constraints.m - Contains a function to compute the equilibrium constraints of the model. Used by KNITRO in optimization and also in solving for equilibrium. Also supplies
the gradient of the constraints. 

eqm_constraints.m -  Contains a function to compute the equilibrium constraints of the model. Used by KNITRO in optimization and also in solving for equilibrium. Also supplies
the gradient of the constraints. The only difference between it and model_constraints is that the parameter vector is not included as part of the input vector. 

likelihood_func.m - Contains a function to compute the model likelihood. Used as objective function by KNITRO in estimation. 

knitro.opt - KNITRO options file. 

INPUT FILES:
german_data_for_matlab_new.out (PROPRIETARY) - Csv containing all the information on each project in Germany, including distances to assembly.  Each line is a single project.
 
danish_data_for_matlab_new.out - Csv containing all the information on each project in Germany, including distances to assembly.  Each line is a single project.

stateDumSorted.csv (PROPRIETARY) - Contains state of origin information for each project. 

dist_matrix_data.mat (PROPRIETARY) - Inter-project distance matrixes for use in the moran test. 

OUTPUT FILES: 

data_structures.mat (PROPRIETARY) - Contains the structures used as inputs to the estimation. 

est_point.mat - Results from the estimation and counterfactuals at the estimated point. Including optimized point and standard errors. 



