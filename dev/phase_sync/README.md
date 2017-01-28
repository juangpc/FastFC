# FastFC
Fast Functional Connectivity Toolbox. This toolbox implements efficient C computing of different FC indices. It does so through C-MEX external interface appi and within Matlab's development environment.

Development branch for Phase Connectivity measures. 

This branch contains the source-code for implementation of PLV, PLI, wPLI and the imaginary part of the Coherence measures. 

All in a single file in order to help ports to other toolboxes. Matlab's MEX external interface has been used to connect Matlab's environment variables to this C function.

There is also a Phase Sync function version which implements FIR filtereing at the same time. However this version is not finished. There is aproblem with it and it blows matlab.

