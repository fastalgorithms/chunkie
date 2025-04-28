%KERNELCLASSTEST check the use of kernel class for sending FMM
% and singularity info to various routines
%
% 

seed = 8675309;
rng(seed);

doadap = false;

% geometry parameters and construction

cparams = [];
cparams.eps = 1.0e-12;
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

% kernel defs

kernd = kernel('lap','d');
kerns = kernel('lap','s');

opdims = [1 1];

% eval u and dudn on boundary

srcinfo = []; srcinfo.r = sources; 

[ubdry,gradubdry] = kerns.fmm(eps,srcinfo,chnkr.r(:,:),strengths);
unbdry = sum(chnkr.n(:,:).*gradubdry,1);

%

% eval u at targets

utarg = kerns.fmm(eps,srcinfo,targets,strengths);


% test green's id

opts = [];
start=tic; Du = chunkerkerneval(chnkr,kernd,ubdry,targets,opts); 
toc(start)
start=tic; Sun = chunkerkerneval(chnkr,kerns,unbdry,targets,opts); 
toc(start)

utarg2 = Sun-Du;

%

relerr = norm(utarg-utarg2,'fro')/norm(utarg,'fro');

fprintf('relative frobenius error %5.2e\n',relerr);

assert(relerr < 1e-11);

nankern = kernel.nans();
assert(nankern.isnan);

kerntmp = kernd + nankern;
assert(kerntmp.isnan);

kerntmp = nan*kerns;
assert(kerntmp.isnan);

