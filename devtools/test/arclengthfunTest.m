%ARCLENGTHFUNTEST tests the arclengthfun routine for computing distance
% along length of curve
% 

seed = 8675309;
rng(seed);

% geometry parameters and construction


cparams = [];
cparams.eps = 1.0e-9;
pref = []; 
pref.k = 16;
narms = 0;
amp = 0.0;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref);

[s, ~, chnkr] = arclengthfun(chnkr);
% Compare to analytic arclength for circle
ts = squeeze(atan2(chnkr.r(2,:,:), chnkr.r(1,:,:)));
ts(ts<0) = ts(ts<0) + 2*pi;
assert(norm(s-ts) < 1e-12);

% Now test two circles
chnkrs(2) = chunker();
chnkrs(1) = chnkr;
rfac = 1.1;
chnkrs(2) = move(chnkr, [0;0], [3;0], 0, rfac);
chnkrtotal = merge(chnkrs);
[s, nchs, ~] = arclengthfun(chnkrtotal);

assert(norm(s(:,1:nchs(1)) - ts) < 1e-12);
assert(norm(s(:,nchs(1)+1:end) - rfac*ts) < 1e-12);
