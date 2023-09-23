%CHUNKERINTERIORTEST tests the routines for checking whether a point is 
% inside a domain or not
% 

clearvars; close all;
seed = 8675309;
rng(seed);
addpaths_loc();

doadap = false;

% geometry parameters and construction


cparams = [];
cparams.eps = 1.0e-9;
pref = []; 
pref.k = 16;
narms = 5;
amp = 0.5;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

% nxdir = 20;
% rmin = min(chnkr);
% rmax = max(chnkr);
% xgrid = linspace(rmin(1),rmax(1),nxdir);
% ygrid = linspace(rmin(2),rmax(2),nxdir);
% [xx,yy] = meshgrid(xgrid,ygrid);
% nt = length(xx(:)); targs = zeros(2,nt); 
% targs(1,:) = xx(:); targs(2,:) = yy(:);

nt = 10000;
scal = 2*rand(1,nt);
tr = 2*pi*rand(1,nt);

targs = bsxfun(@times,starfish(tr,narms,amp),scal);


opts = [];
opts.flam = false;
opts.fmm = false;
start = tic; in = chunkerinterior(chnkr,targs,opts); t1 = toc(start);

opts = [];
opts.fmm = false;
opts.flam = true;
start = tic; in2 = chunkerinterior(chnkr,targs,opts); t2 = toc(start);

opts = [];
opts.fmm = true;
opts.flam = false;
start = tic; in3 = chunkerinterior(chnkr,targs,opts); t3 = toc(start);


fprintf('%5.2e s : time for chunkerinterior (no flam)\n',t1);
fprintf('%5.2e s : time for chunkerinterior (with flam)\n',t2);
fprintf('%5.2e s : time for chunkerinterior (with fmm)\n',t3);

assert(all(in(:) == (scal(:) < 1)));
assert(all(in2(:) == (scal(:) < 1)));
assert(all(in3(:) == (scal(:) < 1)));

x1 = linspace(-2,2,100);
[xx,yy] = meshgrid(x1,x1);

opts = [];
opts.fmm = true;
opts.flam = false;
start = tic; in3 = chunkerinterior(chnkr,{x1,x1},opts); t4 = toc(start);
fprintf('%5.2e s : time for chunkerinterior (meshgrid, with fmm)\n',t4);