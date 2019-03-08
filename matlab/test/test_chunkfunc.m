
%TEST_CHUNKFUNC
%
% This file tests the routine chunkfunc on a couple of examples

addpath('../src','../../mwrap')

cparams.eps = 1.0e-6;
cparams.nchmax = 100000;
narms = 10;
amp = 0.5;
start = tic; chunker = chunkfunc(@(t) starfish(t,narms,amp),cparams); 
toc(start)

xs = chunker.chunks(1,:,:); xs = xs(:);
ys = chunker.chunks(2,:,:); ys = ys(:);

figure(1)
clf
subplot(1,2,1)
scatter(xs,ys)
axis equal

modes = randn(11,1); modes(1) = 1.1*sum(abs(modes(2:end))); ctr = [1.0;-0.5];

start = tic; chunker = chunkfunc(@(t) curvebymode(t,modes,ctr),cparams); 
toc(start)

xs = chunker.chunks(1,:,:); xs = xs(:);
ys = chunker.chunks(2,:,:); ys = ys(:);

subplot(1,2,2)
scatter(xs,ys)
axis equal
