chunkerkerneval_greenlapTest0();


function chunkerkerneval_greenlapTest0()
%CHUNKERKERNEVAL_GREENLAPTEST test the routines for integrating over 
% chunks against the Green's ID for Laplace
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

ns = 100;
ts = 0.0+2*pi*rand(ns,1);
sources = starfish(ts,narms,amp);
sources = 3.0*sources;
strengths = randn(ns,1);

% targets

nt = 300;
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
kernsprime = kernel('lap','sprime');

opdims = [1 1];

% eval u and dudn on boundary

srcinfo = []; srcinfo.r = sources; 
targinfo = []; targinfo.r = chnkr.r(:,:); 
targinfo.d = chnkr.d(:,:);
kernmats = kerns.eval(srcinfo,targinfo);
kernmatsprime = kernsprime.eval(srcinfo,targinfo);
densu = kernmats*strengths;
densun = kernmatsprime*strengths;

% eval u at targets

targinfo = []; targinfo.r = targets;
kernmatstarg = kerns.eval(srcinfo,targinfo);
utarg = kernmatstarg*strengths;


% test green's id. we use FLAM, direct and FMM for smooth work

opts = [];
opts.flam = true;
start=tic; Du = chunkerkerneval(chnkr,kernd,densu,targets,opts); 
toc(start)
start=tic; Sun = chunkerkerneval(chnkr,kerns,densun,targets,opts); 
toc(start)

utarg2 = Sun-Du;

%

relerr = norm(utarg-utarg2,'fro')/norm(utarg,'fro');

fprintf('relative frobenius error, forcing flam %5.2e\n',relerr);

assert(relerr < 1e-11);

opts = [];
opts.accel = false;
start=tic; Du = chunkerkerneval(chnkr,kernd,densu,targets,opts); 
toc(start)
start=tic; Sun = chunkerkerneval(chnkr,kerns,densun,targets,opts); 
toc(start)

utarg2 = Sun-Du;

%

relerr = norm(utarg-utarg2,'fro')/norm(utarg,'fro');

fprintf('relative frobenius error, forcing direct %5.2e\n',relerr);

assert(relerr < 1e-11);

opts = [];
opts.forcefmm = true;
start=tic; Du = chunkerkerneval(chnkr,kernd,densu,targets,opts); 
toc(start)
start=tic; Sun = chunkerkerneval(chnkr,kerns,densun,targets,opts); 
toc(start)

utarg2 = Sun-Du;

%

relerr = norm(utarg-utarg2,'fro')/norm(utarg,'fro');

fprintf('relative frobenius error, forcing fmm %5.2e\n',relerr);

assert(relerr < 1e-11);



end


