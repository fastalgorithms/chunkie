
%TESTCHUNKSLIB
%
% This file tests whether or not chunks 
% lib is wrapped properly

cparams.eps = 1.0e-6;
cparams.nchmax = 100000;
narms = 10;
amp = 0.5;
start = tic; chunker = chunkfunc(@(t) starfish(t,narms,amp),cparams); toc(start)

chunker.nch

xs = chunker.chunks(1,:,:); xs = xs(:);
ys = chunker.chunks(2,:,:); ys = ys(:);

figure(1)
clf
scatter(xs,ys)
axis equal
