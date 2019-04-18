
%TEST_CHUNKFUNC
%
% This file tests the routine chunkfunc on a couple of examples
% as well as testing the plot, quiver, sort, and reverse utilities

addpaths_loc();
cparams = [];
cparams.eps = 1.0e-4;
pref = []; 
pref.k = 16;
narms = 10;
amp = 0.5;
start = tic; chnkr = chunkfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

%

fprintf('%5.2e seconds to chunk starfish with %d chunks\n',t1,chnkr.nch);

figure(1)
clf
subplot(1,2,1)
plot(chnkr)
hold on
quiver(chnkr)
axis equal

%

modes = randn(11,1); modes(1) = 1.1*sum(abs(modes(2:end))); ctr = [1.0;-0.5];

start = tic; chnkr = chunkfunc(@(t) curvebymode(t,modes,ctr),cparams); 
t1 = toc(start);

fprintf('%5.2e seconds to chunk random mode domain with %d chunks\n', ...
    t1,chnkr.nch);

chnkr = reverse(chnkr);

subplot(1,2,2)
plot(chnkr)
hold on
quiver(chnkr)
axis equal
