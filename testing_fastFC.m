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

y_filtfilt=filtfilt(b,1,x);
tic;y_filtfilt=filtfilt(b,1,x);t_filtfilt=toc;
y_fastfc=fastfc_filt(b,x,1); %this first call takes longer for optimization purposes.
tic;y_fastfc=fastfc_filt(b,x,1);t_fastfc=toc;
disp(' New filtfilt version called ''fastfc_filt'' is working properly.');
disp([' FastFC is ' num2str(t_filtfilt/t_fastfc) ' times faster than Matlab''s filtfilt.']);
disp(' ');

%% Testing Phase Synchronization function
disp(' Starting Phase Synchronization test:');

tic;[plv,pval_plv,pli,wpli]=fastfc_ps(x,n_samples*.1,1);t_ps=toc;

disp(' Function ''fastfc_ps'' is working properly.');
disp([' FastFC calculation of Phase Synch. indices: PLV, PLI and wPLI took: ' num2str(t_ps) ' seconds.']);
disp(' ');

%% Testing Mutual Info
disp(' Starting Mutual Information test:');

emb_dim=3;
tau=2;
k=4;
tic;[mi]=fastfc_mi(x,emb_dim,tau,k);t_mi=toc;

disp(' Function ''fastfc_mi'' is working properly.');
disp([' FastFC calculation of Mutual Information index took: ' num2str(t_mi) ' seconds.']);
disp(' ');

%% Testing Generalized Synchronization
disp(' Starting Generalized Synchronization test:');

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
