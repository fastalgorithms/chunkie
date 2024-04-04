%CHUNKERINTEGRALTEST test the routines for integrating over chunks
%
% 

clearvars; clear all;
addpaths_loc();

seed = 8675309;
rng(seed);

% geometry parameters and construction

cparams = [];
cparams.eps = 1.0e-4;
pref = []; 
pref.k = 16;
narms = 5;
amp = 0.5;
chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams);


% build helmholtz dirichlet matrix

fkern = @(s,t) chnk.helm2d.kern(zk,s,t,'D');
