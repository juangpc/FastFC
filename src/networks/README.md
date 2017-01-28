# FastFC
Fast Functional Connectivity Toolbox. This toolbox implements efficient C computing of different FC indices. It does so through C-MEX external interface appi and within Matlab's development environment.

This branch is dedicated to the development of network indices, graph theoretic measures of FC adjacency matrices.
Be aware autolinks are forbidden. Thus, adjacency matrices' principal diagonal must be zeroed. 
## Contents
- Clustering Coefficient for binary undirected networks. 
- Clustering Coefficient for weighted undirected networks. 
- Strength for weighted undirected networks. 
- Shortest path length for weighted networks (both directed and undirected). 
- Betweenness centrality for weighted networks (both directed and undirected). 

## Future Development
- Shortest path length could perform much better if the inversion of the adj. matrix could be done by thread parallelization. At the same time double to float transformation could take place. And continue there as the program does now.

