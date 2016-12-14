# FastFC
Fast Functional Connectivity Toolbox. This toolbox implements efficient C computing of different FC indices. It does so through C-MEX external interface appi and within Matlab's development environment.

This is the main release branch of the repo. Please see other branches for accesing source code. 

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



