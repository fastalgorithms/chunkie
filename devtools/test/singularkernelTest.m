%SINGULARKERNELTEST check that principal value and hypersingular type
% quadratures are working.
% 

clearvars; close all;
seed = 8675309;
rng(seed);
addpaths_loc();

doadap = false;

% geometry parameters and construction

cparams = [];
cparams.eps = 1.0e-9;
pref = []; 
pref.k = 16;
narms = 5;
amp = 0.3;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

% sources

ns = 10;
ts = 0.0+2*pi*rand(ns,1);
sources = starfish(ts,narms,amp);
sources = 3.0*sources;
strengths = randn(ns,1);

hold off
plot(chnkr)
hold on
scatter(sources(1,:),sources(2,:),'o')
axis equal 

%

% kernel defs

kernd = kernel('lap','d');
kerns = kernel('lap','s');
kerndprime = kernel('lap','dprime');
kernsprime = kernel('lap','sprime');
kernstau = kernel('lap','stau');

% eval u and dudn on boundary

srcinfo = []; srcinfo.r = sources; 

eps = 1e-15;
[ubdry,gradubdry] = kerns.fmm(eps,srcinfo,chnkr.r(:,:),strengths,2);
unbdry = sum(chnkr.n(:,:).*gradubdry,1);
tau = -chnk.perp(chnkr.n(:,:));
utbdry = sum(tau(:,:).*gradubdry,1);

wtsc = weights(chnkr);
sprimemat = chunkermat(chnkr,kernsprime);
staumat = chunkermat(chnkr,kernstau);
sys = 0.5*eye(chnkr.npt) + sprimemat + ones(chnkr.npt,1)*(wtsc(:).');
mu = sys\unbdry(:);
utau = staumat*mu;

relerr = norm(utau-utbdry(:))/norm(utbdry);
fprintf('%5.2e : error in utau test (testing pv integration)\n',relerr);
assert(relerr < 1e-9);

dmat = chunkermat(chnkr,kernd);
dprimemat = chunkermat(chnkr,kerndprime);
sys = -0.5*eye(chnkr.npt) + dmat;
mu = sys\ubdry(:);
un = dprimemat*mu;

relerr = norm(un-unbdry(:))/norm(unbdry);
fprintf('%5.2e : error in un test (testing hs integration)\n',relerr);
