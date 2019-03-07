%TEST_GREENLAP test the routines for integrating over chunks against the
% Green's id for Laplace
%
% 

seed = 8675309;
rng(seed);
addpath('../src')
addpath('../../mwrap')

doadap = false;

% geometry parameters and construction

cparams.eps = 1.0e-6;
cparams.nchmax = 100000;
cparams.nover = 1;
narms = 5;
amp = 0.5;
chunker = chunkfunc(@(t) starfish(t,narms,amp),cparams);

% sources

ns = 10;
ts = 0.0+2*pi*rand(ns,1);
fvals = starfish(ts,narms,amp);
sources = zeros(2,ns);
sources(1,:) = fvals(:,1); sources(2,:) = fvals(:,2);
sources = 1.5*sources;
strengths = randn(ns,1);

% targets

nt = 100;
ts = 0.0+2*pi*rand(nt,1);
fvals = starfish(ts,narms,amp);
targets = zeros(2,nt);
targets(1,:) = fvals(:,1); targets(2,:) = fvals(:,2);
targets = targets.*repmat(rand(1,nt),2,1);

% plot geo and sources

xs = chunker.chunks(1,:,:); xmin = min(xs(:)); xmax = max(xs(:));
ys = chunker.chunks(2,:,:); ymin = min(ys(:)); ymax = max(ys(:));

hold off
scatter(xs(:),ys(:))
hold on
scatter(sources(1,:),sources(2,:),'o')
scatter(targets(1,:),targets(2,:),'x')

% kernel defs

kernd = @(s,t,sn,tn) glapkern(s,t,sn,tn,'d');
kerns = @(s,t,sn,tn) glapkern(s,t,sn,tn,'s');
kernsprime = @(s,t,sn,tn) glapkern(s,t,sn,tn,'sprime');

% eval u and dudn on boundary

targs = chunker.chunks; targs = reshape(targs,2,chunker.k*chunker.nch);
targsn = chunknormals(chunker); 
targsn = reshape(targsn,2,chunker.k*chunker.nch);

kernmats = kerns(sources,targs,[],targsn);
kernmatsprime = kernsprime(sources,targs,[],targsn);
densu = kernmats*strengths;
densun = kernmatsprime*strengths;

% eval u at targets

kernmatstarg = kerns(sources,targets,[],[]);
utarg = kernmatstarg*strengths;


% test green's id

opts.usesmooth=false;
opts.verb=false;
opts.quadkgparams = {'RelTol',1.0e-8,'AbsTol',1.0e-8};
start=tic; Du = chunkerintkern(chunker,kernd,ndims,densu,targets,opts); 
toc(start)
start=tic; Sun = chunkerintkern(chunker,kerns,ndims,densun,targets,opts); 
toc(start)

utarg2 = Sun-Du;

%

norm(utarg-utarg2,'fro')

