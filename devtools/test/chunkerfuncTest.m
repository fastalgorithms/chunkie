
%CHUNKERFUNCTEST
%
% This file tests the routine chunkerfunc on a couple of examples

%

clearvars; close all;
addpaths_loc();
cparams = [];
cparams.eps = 1.0e-4;
pref = []; 
pref.k = 16;
narms = 10;
amp = 0.5;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

[~,~,info] = sortinfo(chnkr);
assert(info.ier == 0,'adjacency issues after chunk build starfish');

%

fprintf('%5.2e seconds to chunk starfish with %d chunks\n',t1,chnkr.nch);

cparams.nout = 3;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

[~,~,info] = sortinfo(chnkr);
assert(info.ier == 0,'adjacency issues after chunk build starfish');

%

fprintf('%5.2e seconds to chunk starfish with %d chunks\n',t1,chnkr.nch);
% 
% figure(1)
% clf
% subplot(1,3,1)
% plot(chnkr)
% hold on
% quiver(chnkr)
% axis equal

% chunk up random curve

modes = randn(11,1); modes(1) = 1.1*sum(abs(modes(2:end))); ctr = [1.0;-0.5];

start = tic; chnkr = chunkerfunc(@(t) chnk.curves.bymode(t,modes,ctr),cparams); 
t1 = toc(start);

fprintf('%5.2e seconds to chunk random mode domain with %d chunks\n', ...
    t1,chnkr.nch);

[~,~,info] = sortinfo(chnkr);
assert(info.ier == 0,'adjacency issues after chunk build by modes');

chnkr = reverse(chnkr);

[~,~,info] = sortinfo(chnkr);
assert(info.ier == 0,'adjacency issues after chunker reversal');


% subplot(1,3,2)
% plot(chnkr)
% hold on
% quiver(chnkr)
% axis equal

% chunk up circle and test area

r = 5*rand(); ctr = [1.0;-0.5];

circfun = @(t) ctr + r*[cos(t(:).');sin(t(:).')];
start = tic; chnkr = chunkerfunc(circfun,cparams); 
t1 = toc(start);

fprintf('%5.2e seconds to chunk circle domain with %d chunks\n', ...
    t1,chnkr.nch);

[~,~,info] = sortinfo(chnkr);
assert(info.ier == 0,'adjacency issues after chunk build circle');

a = area(chnkr);
assert(abs(a - pi*r^2) < 1e-12,'area wrong for circle domain')


% subplot(1,3,3)
% plot(chnkr)
% hold on
% quiver(chnkr)
% axis equal
