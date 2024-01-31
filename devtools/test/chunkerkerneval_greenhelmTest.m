%CHUNKERKERNEVAL_GREENHELMTEST test the routines for integrating over 
% chunks against the Green's ID for Helmholtz
%
% 

clearvars; close all;
seed = 8675309;
rng(seed);
addpaths_loc();

% geometry parameters and construction

cparams = [];
cparams.eps = 1.0e-11;
pref = []; 
pref.k = 16;
narms = 5;
amp = 0.5;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

% sources

ns = 10;
ts = 0.0+2*pi*rand(ns,1);
sources = starfish(ts,narms,amp);
sources = 3.0*sources;
strengths = randn(ns,1);

% targets

nt = 100;
ts = 0.0+2*pi*rand(nt,1);
targets = starfish(ts,narms,amp);
targets = targets.*repmat(rand(1,nt),2,1);

% plot geo and sources

xs = chnkr.r(1,:,:); xmin = min(xs(:)); xmax = max(xs(:));
ys = chnkr.r(2,:,:); ymin = min(ys(:)); ymax = max(ys(:));

hold off
plot(chnkr)
hold on
scatter(sources(1,:),sources(2,:),'o')
scatter(targets(1,:),targets(2,:),'x')
axis equal 

%

zk = rand() + 1i*rand();

% kernel defs

kernd = kernel('h','d',zk);
kerns = kernel('h','s',zk);
kernsprime = kernel('h','sprime',zk);

opdims = [1 1];

% eval u and dudn on boundary

srcinfo = []; srcinfo.r = sources;
targinfo = []; targinfo.r = chnkr.r(:,:); 
targinfo.d = chnkr.d(:,:);
targinfo.n = chnkr.n(:,:);
[u,gradu] = kerns.fmm(1e-12,srcinfo,targinfo,strengths);
densu = u;
densun = sum(chnkr.n(:,:).*gradu,1);

% eval u at targets

targinfo = []; targinfo.r = targets;
utarg = kerns.fmm(1e-12,srcinfo,targinfo,strengths);


% test green's id

start=tic; Du = chunkerkerneval(chnkr,kernd,densu,targets); 
toc(start)
start=tic; Sun = chunkerkerneval(chnkr,kerns,densun,targets); 
toc(start)

utarg2 = Sun-Du;

%

relerr = norm(utarg-utarg2,'fro')/norm(utarg,'fro');

fprintf('relative frobenius error %5.2e\n',relerr);

assert(relerr < 1e-11);

