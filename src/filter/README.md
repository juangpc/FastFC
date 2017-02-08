# FastFC
Fast Functional Connectivity Toolbox. This toolbox implements efficient C computing of different FC indices. It does so through C-MEX external interface appi and within Matlab's development environment.

Development branch for Filtering with zero phase distortion. 

This branch contains the source-code for implementation of a zero phase distortion filtering procedure. 

All in a single file in order to help ports to other toolboxes. Matlab's MEX external interface has been used to connect Matlab's environment variables to this C function.

fastFC_filt function is developed after matlab's filtfilt function. It implements a two passes FIR filtering with previous circular padding of the original signal.

# Related research:
Efficient computation of functional brain networks: towards real-time functional connectivity
Frontiers in Neuroinformatics (2017) García-Prieto Juan, Bajo Ricardo, Pereda Ernesto
