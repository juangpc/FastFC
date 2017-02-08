# FastFC
Fast Functional Connectivity Toolbox
#
This branch is only for developing of Mutual Information index.

This implementation is built on top of MILCA implementation 
(see: http://www.ucl.ac.uk/ion/departments/sobell/Research/RLemon/MILCA/MILCA)

Here you will find: 

- A modified version of binutils in order to be compatible with windows Visual Studio compatible version of binutils.cpp and binutils.h
- A binary verson of MILCA but compiled directly with Visual Studio.
- A Matlab version of Mutual Information index calculation which calls the binary program just preiously mentioned.
- The source coude of a modified version of MILCA program but implemented so that it communicate with MATLAB through its external MEX interface, instead of through a text document.
- A binary version of this last program, compiled with Visual Studio for Matlab in Windows 64bit systems.
- Some other tools and tries that have been generated while development. See log for more details.

# Related research:
Efficient computation of functional brain networks: towards real-time functional connectivity
Frontiers in Neuroinformatics (2017) García-Prieto Juan, Bajo Ricardo, Pereda Ernesto
