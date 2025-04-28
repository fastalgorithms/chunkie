
%CHNKFUNCUNITEST
%
% This file tests the routine chunkerfunc on a couple of examples
% as well as testing the plot, quiver, sort, and reverse utilities


nch = 10;

narms = 3;
amp = 0.5;
start = tic; chnkr = chunkerfuncuni(@(t) starfish(t,narms,amp),nch); 
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

start = tic; chnkr = chunkerfuncuni(@(t) chnk.curves.bymode(t,modes,ctr),nch);
t1 = toc(start);

fprintf('%5.2e seconds to chunk random mode domain with %d chunks\n', ...
    t1,chnkr.nch);

chnkr = reverse(chnkr);

subplot(1,2,2)
plot(chnkr)
hold on
quiver(chnkr)
axis equal


% chunk up circle and test area

r = 5*rand(); ctr = [1.0;-0.5];

circfun = @(t) ctr + r*[cos(t(:).');sin(t(:).')];

start = tic; chnkr = chunkerfuncuni(circfun,nch); 
t1 = toc(start);

fprintf('%5.2e seconds to chunk circle domain with %d chunks\n', ...
    t1,chnkr.nch);

[~,~,info] = sortinfo(chnkr);
assert(info.ier == 0,'adjacency issues after chunk build circle');

a = area(chnkr);
assert(abs(a - pi*r^2) < 1e-12,'area wrong for circle domain')


