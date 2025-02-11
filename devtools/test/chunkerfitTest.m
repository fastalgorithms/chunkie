
%CHUNKERFITTEST

clearvars; close all;
addpaths_loc();

% Sample a smooth curve at random points
rng(0)
n = 20;
tt = sort(2*pi*rand(n,1));
r = chnk.curves.bymode(tt, [2 0.5 0.2 0.7]);

opts = [];
opts.ifclosed = true;
opts.cparams = [];
opts.cparams.eps = 1e-6;
opts.pref = [];
opts.pref.k = 16;
chnkr = chunkerfit(r, opts);
assert(checkadjinfo(chnkr) == 0);

opts.ifclosed = false;
chnkr = chunkerfit(r(:,1:10), opts);
assert(checkadjinfo(chnkr) == 0);
