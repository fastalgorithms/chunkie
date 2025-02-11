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


% scalar function on boundary

fscal = @(xx) cos(xx(1,:)-1.0) + sin(xx(2,:)-0.5);
fvals = fscal(reshape(chnkr.r,2,chnkr.k*chnkr.nch));
opts = [];
opts.quadgkparams = {'RelTol',1e-15};
opts.usesmooth = false;
fscal_int1 = chunkerintegral(chnkr,fvals,opts);
opts.usesmooth = true;
fscal_int3 = chunkerintegral(chnkr,fvals,opts);

opts.usesmooth = false;
fscal_int2 = chunkerintegral(chnkr,fscal,opts);
opts.usesmooth = true;
fscal_int4 = chunkerintegral(chnkr,fscal,opts);

assert(abs(fscal_int1-fscal_int2)/abs(fscal_int2) < 1e-9);
assert(abs(fscal_int3-fscal_int2)/abs(fscal_int2) < 1e-9);
assert(abs(fscal_int4-fscal_int2)/abs(fscal_int2) < 1e-9);