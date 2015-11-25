%    This file is part of Fast Functional Connectivity (FastFC)
%
%    FastFC is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FastFC is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FastFC.  If not, see <http://www.gnu.org/licenses/>.
%
%    ------------------------------------------ 
%    Contact: Juan Garcia-Prieto    juangpc (at) gmail.com
%    ------------------------------------------
%
%    Please consider helping by citing our research.
% 
%    J. Garcia-Prieto, E. Pereda
%

%% FastFC testing script

% configure number of sensors and samples
n_samples=1000;
n_sensors=100;
x=randn(n_samples,n_sensors); % remember columnwise matrix of each sensor.

%% 
disp(' '); 
disp(' '); 
disp(' FastFC testing script');
disp(' ');
disp(' ');
disp(' This function is intended to check the correct installation of FastFC. ');
disp(' While it helps understanding how different functions can be used. ');
disp(' ');
disp([' The testing setup will be a set of ' num2str(n_samples) ' samples and ' num2str(n_sensors) ' sensors.']);
disp(' ..................................');

disp(' ');
% 
%% testing FastFC filtering function
disp(' Starting Zero Phase Distortion filtering test:')
b=fir1(ceil(n_samples/5),[0.0160 0.0240],'bandpass',hamming(ceil(n_samples/5)+1),'scale');

% y=fastfc_filt(filter,data,mode)
% Input parameters:
% filter = row vector with filter denominator coefficients.
% data = column matrix with data to be filtered. Each sensor should be each column.
% mode = mode for FFT scheduling
% mode = 0 ,  execute fastest but suboptimal.
% mode = 1 , execute fast but in a suboptimal algorithm.
% mode = 2 , execute slower the first time, but consider possible faster algorithms.
% mode = 3 , execute slowest the first time, but consider the fastest algorithm.
% Output parameters:
% y = array the same size of data, with each column corresponding to each filtered column of data.

y_filtfilt=filtfilt(b,1,x);
tic;y_filtfilt=filtfilt(b,1,x);t_filtfilt=toc;
y_fastfc=fastfc_filt(b,x,1); %this first call takes longer for optimization purposes.
tic;y_fastfc=fastfc_filt(b,x,1);t_fastfc=toc;
disp(' New filtfilt version called ''fastfc_filt'' is working properly.');
disp([' FastFC is ' num2str(t_filtfilt/t_fastfc) ' times faster than Matlab''s filtfilt.']);
disp(' ');

%% Testing Phase Synchronization function
disp(' Starting Phase Synchronization test:');

% [plv,pval_plv,pli,wpli]=fastfc_ps(data,samples_to_discard,mode)
% 
% Input parameters:
% data = data (sensors by columns).
% samples_to_discard = samples to discard at the beginning and samples to discard at the end of the phase signals. i.e. To discard 200 samples at the beginning and 200 samples at the end of each sensor, samples_to_discard should equal 200.
%         mode = 0 -> execute fastest but suboptimal.
%         mode = 1 -> execute fast but in a suboptimal algorithm.
%         mode = 2 -> execute slower the first time, but consider possible faster algorithms.
%         mode = 3 -> execute slowest the first time, but consider the fastest algorithm.
% Output parameters:
% plv = Phase Locking Value Functional Connectivity matrix.
% pval_plv = pvalue for each index of the PLV matrix. 
% pli = Phase Locking Index Functional Connectivity matrix.
% wPli = weighted Phase Locking Functional Connectivity matrix.

tic;[plv,pval_plv,pli,wpli]=fastfc_ps(x,n_samples*.1,1);t_ps=toc;

disp(' Function ''fastfc_ps'' is working properly.');
disp([' FastFC calculation of Phase Synch. indices: PLV, PLI and wPLI took: ' num2str(t_ps) ' seconds.']);
disp(' ');

%% Testing Mutual Info
disp(' Starting Mutual Information test:');

% [mi]=fastfc_mi(data,emb_dim,tau,k);
% Input parameters:
% data = eeg data (sensors by columns).
% emb_dim = embedding dimensions to consider
% tau = time lag to consider for embedding
% k = number of neighbours to consider
% Output parameters:
% mi = Mutual Information Functional Connectivity matrix
	
emb_dim=3;
tau=2;
k=4;
tic;[mi]=fastfc_mi(x,emb_dim,tau,k);t_mi=toc;

disp(' Function ''fastfc_mi'' is working properly.');
disp([' FastFC calculation of Mutual Information index took: ' num2str(t_mi) ' seconds.']);
disp(' ');

%% Testing Generalized Synchronization
disp(' Starting Generalized Synchronization test:');

% [S,H,M,L]=fastfc_gs(data,emb_dim,tau,k,w,states_eff_step)
% 
% Input parameters:
% data = eeg data (sensors by columns).
% emb_dim = embedding dimension 
% tau = time lag for embedding
% k = number of neighbours to consider
% w = window correction for neighbour finding
% states_eff_step = state-space down sampling to consider when calculating distances 
% Output parameters: 
% S, H, M and L = Functional Connectivity matrices for each index.

emb_dim=3;
tau=2;
k=4;
w=20;
states_eff_step=1;

tic;[S,H,M,L]=fastfc_mi(x,emb_dim,tau,k,w,states_eff_step);t_gs=toc;

disp(' Function ''fastfc_gs'' is working properly.');
disp([' FastFC calculation of Generalized Synchronization indices: S, H, M, and L took: ' num2str(t_gs) ' seconds.']);
disp(' ');

%% Testing Strength
disp(' Starting Strength test:');

A=randn(n_sensors);
A(1:n_sensors+1:end)=0;

% [S]=fastfc_strength_wu(A)
% 
% Input parameters: 
% A = adjacency matrix of real values between 0 and 1, with zeroed principal diagonal elements.
% Output parameters:
% C = a row matrix where every value represents the Strength of each node. 

tic;[S]=fastfc_strength_wu(A);t_S=toc;

disp(' Function ''fastfc_strength_wu'' is working properly.');
disp([' FastFC calculation of Strength index took: ' num2str(t_S) ' seconds.']);
disp(' ');

%% Testing Clustering Coefficient
disp(' Starting Clustering Coefficient test:');

A=randn(n_sensors);
A(1:n_sensors+1:end)=0;

tic;[C]=fastfc_cluster_coef_wu(A);t_C=toc;

disp(' Function ''fastfc_cluster_coef_wu'' is working properly.');
disp([' FastFC calculation of Clustering Coefficient index took: ' num2str(t_C) ' seconds.']);
disp(' ');

%% Testing Shortest Path Length
disp(' Starting Shortest Path Lenght test:');

A=randn(n_sensors);
A(1:n_sensors+1:end)=0;
W=(A.^-1);

tic;[D,L]=fastfc_shortest_path_length_w(W);t_L=toc;

disp(' Function ''fastfc_shortest_path_length_w'' is working properly.');
disp([' FastFC calculation of Shortest Path Length index took: ' num2str(t_L) ' seconds.']);
disp(' ');

%% Testing Betweenness Centrality
disp(' Starting Betweenness Centrality test:');

A=randn(n_sensors);
A(1:n_sensors+1:end)=0;
W=(A.^-1);

tic;[D,L,B]=fastfc_betweenness_cent_w(W);t_B=toc;

disp(' Function ''fastfc_shortest_path_length_w'' is working properly.');
disp([' FastFC calculation of Shortest Path Length index took: ' num2str(t_B) ' seconds.']);
disp(' ');

%% Final
disp(' ............ Fast FC is working rock solid. Thanks! ..............')
disp(' ');
disp(' This file is part of Fast Functional Connectivity (FastFC)');
disp(' Please consider helping by citing our work');
disp(' http://juangpc.github.io/FastFC/');
disp(' ');
