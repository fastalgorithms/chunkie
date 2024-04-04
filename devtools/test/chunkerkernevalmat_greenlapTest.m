%CHUNKERKERNEVAL_GREENLAPTEST test the routines for integrating over 
% chunks against the Green's ID for Laplace
%
% 

clearvars; close all;
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

kernd = @(s,t) chnk.lap2d.kern(s,t,'d');
kerns = @(s,t) chnk.lap2d.kern(s,t,'s');
kernsprime = @(s,t) chnk.lap2d.kern(s,t,'sprime');

opdims = [1 1];

% eval u and dudn on boundary

srcinfo = []; srcinfo.r = sources; 
targinfo = []; targinfo.r = chnkr.r(:,:); 
targinfo.d = chnkr.d(:,:);
kernmats = kerns(srcinfo,targinfo);
kernmatsprime = kernsprime(srcinfo,targinfo);
densu = kernmats*strengths;
densun = kernmatsprime*strengths;

% eval u at targets

targinfo = []; targinfo.r = targets;
kernmatstarg = kerns(srcinfo,targinfo);
utarg = kernmatstarg*strengths;


% test green's id

opts=[];
start=tic; Dmat = chunkerkernevalmat(chnkr,kernd,targets,opts); 
Du = Dmat*densu;
toc(start)
opts=[]; opts.forceadap=true;
start=tic; Dmat = chunkerkernevalmat(chnkr,kernd,targets,opts); 
Du2 = Dmat*densu;
toc(start)
opts=[];
start=tic; Smat = chunkerkernevalmat(chnkr,kerns,targets,opts); 
Sun = Smat*densun;
toc(start)

utarg2 = Sun-Du;
utarg3 = Sun-Du2;

%

relerr = norm(utarg-utarg2,'fro')/norm(utarg,'fro');
assert(relerr < 1e-11);

fprintf('relative frobenius error %5.2e\n',relerr);

relerr = norm(utarg-utarg3,'fro')/norm(utarg,'fro');
assert(relerr < 1e-11);


