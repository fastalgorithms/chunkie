%TEST_INTROUTS test the routines for integrating over chunks
%
% 

seed = 8675309;
rng(seed);
addpath('../src')
addpath('../../mwrap')

% geometry parameters and construction

cparams.eps = 1.0e-6;
cparams.nchmax = 100000;
cparams.nover = 1;
narms = 5;
amp = 0.5;
chunker = chunkfunc(@(t) starfish(t,narms,amp),cparams);

% scalar function on boundary

fscal = @(xx) cos(xx(1,:)-1.0) + sin(xx(2,:)-0.5);
fvals = fscal(reshape(chunker.chunks,2,chunker.k*chunker.nch));
fscal_int1 = chunkerint(chunker,fvals);
opts.quadgkparams = {'RelTol',1e-15};
fscal_int2 = chunkerint(chunker,fscal,opts);

abs(fscal_int1-fscal_int2)/abs(fscal_int2)
