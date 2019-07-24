%TEST_CHUNKERINFLAM
%
% determine points inside/outside domain reasonably fast 

iseed = 8675309;
rng(iseed);

addpaths_loc();

cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 1;
pref = []; 
pref.k = 16;
narms = 3;
amp = 0.25;
start = tic; chnkr = chunkfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

wts = whts(chnkr);

fprintf('%5.2e s : time to build geo\n',t1)

% sources

ns = 10;
ts = 0.0+2*pi*rand(ns,1);
sources = starfish(ts,narms,amp);
sources = 3.0*sources;
strengths = randn(ns,1);

% targets

nt = 10000;
ts = 0.0+2*pi*rand(nt,1);
targets = starfish(ts,narms,amp);
stretch = 2*rand(1,nt);
targets = targets.*repmat(stretch,2,1);

% plot geo and sources

figure(1)
clf
hold on
scatter(sources(1,:),sources(2,:),'o')
scatter(targets(1,:),targets(2,:),'x')
plot(chnkr)
axis equal 


%%

targin_true = stretch <= 1;

start = tic; targin = chunkerinflam(chnkr,targets); 
t1 = toc(start);

fprintf('%5.2e s : time to determine if in/out \n',t1)

nfail = nnz(targin(:) ~= targin_true(:));

fprintf('number of targets:      %7d\n',nt)
fprintf('number of bdry pts:     %7d\n',chnkr.npt)
fprintf('number of failures:     %7d\n',nfail)