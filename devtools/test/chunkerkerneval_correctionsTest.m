
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

ns = 10;
ts = 0.0+2*pi*rand(ns,1);
sources = starfish(ts,narms,amp);
sources = 3.0*sources;
strengths = randn(ns,1);

% targets

nt = 20;
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

zk = rand() + 1i*rand();

skern = kernel('helmholtz','s',zk);
dkern = kernel('helmholtz','d',zk);

rhs = skern.eval(struct("r",sources),chnkr)*strengths;
sys = chunkermat(chnkr,dkern);
sys = -0.5*eye(chnkr.npt)+sys;

sol = sys\rhs;

utrue = skern.eval(struct("r",sources),struct("r",targets))*strengths;

opts = [];
opts.corrections = true;
cormat = chunkerkernevalmat(chnkr,dkern,targets,opts);
opts = [];
opts.cormat = cormat;
u_eval_cor = chunkerkerneval(chnkr,dkern,sol,targets,opts);

assert(norm(utrue-u_eval_cor,inf)<1e-11)


