# FastFC

Fast Functional Connectivity Toolbox. This toolbox implements efficient C computing of different FC indices. It does so through C-MEX external interface appi and within Matlab's development environment. The focus of FastFC is to be very small and easy to adapt to any existing toolbox within Matlab.

This is the main release branch of the repo. It only includes the final executable files, a needed dll and an example script. Please see other branches for accesing source code for each set of functions. 

'Hope it helps!

# Related research:
Efficient computation of functional brain networks: towards real-time functional connectivity
Frontiers in Neuroinformatics (2017) García-Prieto Juan, Bajo Ricardo, Pereda Ernesto

# Installation
The library has been implemented in standard C, although it uses OMP and FFTW (both open source projects). The only thing needed in order to use the functions is the FFTW dll which is included with the executables.

# Filtering tool
- Zero-phase distorsion FIR filtering. 
This filtering function has been built out of necessity.  It closely resembles Matlab's filtfilt built-in function implementation. It performs an N-dimensional circular padding (N=filter order). 

# Functional Connectivity measures
- Phase Synchronization measures: PLV, PLV's p-value, PLI, wPLI and imaginary part of the Coherence.
- Mutual Information. (MILCA implementation).
- Generalized Synchronization measures: S, H, M, L
- Generalized Synchronization measrues (version 2. faster): S, H, M.

# Complex Networks measures
- Strength for weighted undirected networks.
- Clustering Coefficient for binary undirected networks.
- Clustering Coefficient for weighted undirected networks.
- Shortest path length for weighted networks (both directed and undirected).
- Betweenness centrality for weighted networks (both directed and undirected). 



