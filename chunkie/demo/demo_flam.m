%DEMO_FLAM
%
% We have wrappers to using the FLAM library for obtaining a fast-direct
% solver for your system matrix. This is useful for larger scale problems.
%
% This demo calls the FLAM matrix builder and do a basic solve on a 
% triangular chunkgraph

clearvars; close all;
iseed = 8675309;
rng(iseed);

% set up curved triangle geometry

nverts = 3; 
verts = exp(1i*2*pi*(0:(nverts-1))/nverts);
verts = [real(verts);imag(verts)];

edgends = [1:nverts; 2:nverts, 1];
amp = 0.1;
frq = 3;
fchnks    = cell(1,size(edgends,2));
for icurve = 1:size(edgends,2)
    fchnks{icurve} = @(t) sinearc(t,amp,frq);
end
cparams = [];
cparams.nover = 4;
[cgrph] = chunkgraph(verts,edgends,fchnks,cparams);

% sources

ns = 10;
ts = 0.0+2*pi*rand(ns,1);
sources = 3.0*[cos(ts)';sin(ts)'];
strengths = randn(ns,1);

% targets

nt = 1;
ts = 0.0+2*pi*rand(nt,1);
targets = 0.2*[cos(ts)'; sin(ts)'];
targets = targets.*repmat(rand(1,nt),2,1);

% plot geo and sources

figure(1)
clf
hold off
plot(cgrph)
hold on
scatter(sources(1,:),sources(2,:),'o')
scatter(targets(1,:),targets(2,:),'x')
axis equal 

%

kerns = kernel('lap','s');

% eval u on bdry

srcinfo = []; srcinfo.r = sources;
targinfo = []; targinfo.r = cgrph.r(:,:); 
targinfo.d = cgrph.d(:,:);
ubdry = kerns.eval(srcinfo,targinfo)*strengths(:);

% eval u at targets

targinfo = []; targinfo.r = targets;
utarg = kerns.eval(srcinfo,targinfo)*strengths(:);

% to use RCIP properly, kernel must be scaled so system is 
% I + kernel, (I with multiple 1!)

%
% build laplace dirichlet matrix * (-2) densely

kernd = -2*kernel('lap','d');

start = tic; D = chunkermat(cgrph,kernd);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

npt = cgrph.npt;
sys = eye(npt) + D;

% construct FLAM "rskelf" factorization

start = tic; F = chunkerflam(cgrph,kernd,1.0); 
t1 = toc(start);

fprintf('%5.2e s : time to FLAM compress\n',t1)

rhs = ubdry; rhs = rhs(:);
start = tic; sol = gmres(sys,rhs,[],1e-14,100); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

rhs = ubdry; rhs = rhs(:);
start = tic; sol2 = rskelf_sv(F,rhs); t1 = toc(start);

fprintf('%5.2e s : time for rskelf_sv \n',t1)

err = norm(sol-sol2,'fro')/norm(sol,'fro');

fprintf('difference between fast-direct and iterative %5.2e\n',err)

% evaluate at targets and compare

% Collapse cgrph into chnkrtotal
start=tic; Dsol = chunkerkerneval(cgrph,kernd,sol,targets); 
t1 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n',t1)

%

relerr = norm(utarg-Dsol,'fro')/(sqrt(cgrph.npt)*norm(utarg,'fro'));
fprintf('relative frobenius error %5.2e\n',relerr);

%%%
% support functions 

function [r,d,d2] = sinearc(t,amp,frq)
xs = t;
ys = amp*sin(frq*t);
xp = ones(size(t));
yp = amp*frq*cos(frq*t);
xpp = zeros(size(t));
ypp = -frq*frq*amp*sin(t);

r = [(xs(:)).'; (ys(:)).'];
d = [(xp(:)).'; (yp(:)).'];
d2 = [(xpp(:)).'; (ypp(:)).'];
end
