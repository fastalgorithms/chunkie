%GUIDE_KERNELS
%
% This script complements the chunkie guide section on kernel objects.
% It shows the most common methods for building and working with kernels.
%

%%%%%%%%%%%%% defining your own kernel
% START CREATE KERNEL
% start with an empty kernel and add info

kern = kernel();

zk = 1.3; % parameter
kern.eval = @(s,t) expkernel(s,t,zk);
kern.opdims = [1 1]; % for a scalar kernel, dims are [1 1]
% END CREATE KERNEL

% START ADDING KERNELS
% start with an empty kernel and add info

kern1 = kernel('lap','d');
kern2 = kernel('lap','s');

% this creates a combined layer kernel
kern = -2*(kern1 + kern2);
% END ADDING KERNELS

% START MATRIX OF KERNELS
% suppose edge 1 has a dirichlet condition and edge 2 has 
% a neumann condition. we can use a double layer potential 
% on edge 1 and a single layer potential on edge 2

kern11 = kernel('lap','d'); % influence of edge 1 on itself
kern21 = kernel('lap','dprime'); % influence of edge 1 on edge 2
kern12 = kernel('lap','s'); % influence of edge 2 on edge 1
kern22 = kernel('lap','sprime'); % influence of edge 2 on itself

kernmat(2,2) = kernel;
kernmat(1,1) = kern11; kernmat(2,1) = kern21;
kernmat(1,2) = kern12; kernmat(2,2) = kern22;
% END MATRIX OF KERNELS

% START INTERLEAVE
% for a transmission problem it is common to have the solution defined 
% as the sum of two layer potentials. in the exterior of the obstacle, 
% there is one wave speed, k_2 and in the interior there is another, k_1.
% One version of the boundary condition is that the potential and its
% normal derivative should be continuous across the boundary, thus the
% difference in these quantities is prescribed in terms of the incident 
% field

% parameters 
zk1 = 1.5;
zk2 = 1.1;
zks = [zk1,zk2];

% a transmission boundary condition set of kernels 
% first column -> unknown for double layer
% second column -> unknown for single layer
% first row -> jump in potential
% second row -> jump in normal derivative 

% using the "difference" kernels 
kern11 = kernel.helm2ddiff('d',zks); % D_{k1} - D_{k2}
kern21 = kernel.helm2ddiff('dprime',zks); 
kern12 = kernel.helm2ddiff('s',zks); 
kern22 = kernel.helm2ddiff('sprime',zks); 

kernmat(2,2) = kernel;
kernmat(1,1) = kern11; kernmat(2,1) = kern21;
kernmat(1,2) = kern12; kernmat(2,2) = kern22;
kern = kernel(kernmat); % interleaves into a 2x2 matrix-valued kernel
% END INTERLEAVE 