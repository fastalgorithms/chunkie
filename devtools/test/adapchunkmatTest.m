
%TEST_ADAPCHUNKERMAT
%
% define geometry and test adaptive matrix builder routine

iseed = 8675309;
rng(iseed);

addpaths_loc();

zk = randn() + 1i*randn();

cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 2;
pref = []; 
pref.k = 16;
narms = 3;
amp = 0.25;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

% sources

ns = 10;
ts = 0.0+2*pi*rand(ns,1);
sources = starfish(ts,narms,amp);
sources = 3.0*sources;
strengths = randn(ns,1);

% targets

nt = 3;
ts = 0.0+2*pi*rand(nt,1);
targets = starfish(ts,narms,amp);
targets = targets.*repmat(rand(1,nt),2,1);

% plot geo and sources

xs = chnkr.r(1,:,:); xmin = min(xs(:)); xmax = max(xs(:));
ys = chnkr.r(2,:,:); ymin = min(ys(:)); ymax = max(ys(:));

figure(1)
clf
hold off
plot(chnkr)
hold on
scatter(sources(1,:),sources(2,:),'o')
scatter(targets(1,:),targets(2,:),'x')
axis equal 

%

% build layer potential matrix with GGQ routine for comparison

fkern = @(s,t,stau,ttau) chnk.lap2d.kern(s,t,stau,ttau,'S');
fkern = @(s,t,stau,ttau) chnk.helm2d.kern(zk,s,t,stau,ttau,'D');

start = tic; mat1 = chunkermat(chnkr,fkern);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix (GGQ nbor)\n',t1)

%

% use adaptive routine to build matrix (self done by ggq, nbor by adaptive)

type = 'log';
opts = []; opts.robust = false;
opdims = [1 1];
start = tic;
mat2 = chnk.quadadap.buildmat(chnkr,fkern,opdims,type,opts);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix (adap nbor)\n',t1);


% compare

assert(norm(mat1-mat2,'fro')/norm(mat1,'fro') < 1e-9);



