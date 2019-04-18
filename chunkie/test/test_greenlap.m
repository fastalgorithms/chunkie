%TEST_GREENLAP test the routines for integrating over chunks against the
% Green's ID for Laplace
%
% 

seed = 8675309;
rng(seed);
addpaths_loc();

doadap = false;

% geometry parameters and construction

cparams = [];
cparams.eps = 1.0e-12;
pref = []; 
pref.k = 16;
narms = 5;
amp = 0.5;
start = tic; chnkr = chunkfunc(@(t) starfish(t,narms,amp),cparams,pref); 
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

%%

% kernel defs

kernd = @(s,t,sn,tn) glapkern(s,t,sn,tn,'d');
kerns = @(s,t,sn,tn) glapkern(s,t,sn,tn,'s');
kernsprime = @(s,t,sn,tn) glapkern(s,t,sn,tn,'sprime');

ndims = [1 1];

% eval u and dudn on boundary

targs = chnkr.r; targs = reshape(targs,2,chnkr.k*chnkr.nch);
targstau = taus(chnkr); 
targstau = reshape(targstau,2,chnkr.k*chnkr.nch);

kernmats = kerns(sources,targs,[],targstau);
kernmatsprime = kernsprime(sources,targs,[],targstau);
densu = kernmats*strengths;
densun = kernmatsprime*strengths;

% eval u at targets

kernmatstarg = kerns(sources,targets,[],[]);
utarg = kernmatstarg*strengths;


% test green's id

opts.usesmooth=false;
opts.verb=false;
opts.quadkgparams = {'RelTol',1.0e-13,'AbsTol',1.0e-13};
start=tic; Du = chunkerintkern(chnkr,kernd,ndims,densu,targets,opts); 
toc(start)
start=tic; Sun = chunkerintkern(chnkr,kerns,ndims,densun,targets,opts); 
toc(start)

utarg2 = Sun-Du;

%

relerr = norm(utarg-utarg2,'fro')/norm(utarg,'fro');

fprintf('relative frobenius error %5.2e\n',relerr);

