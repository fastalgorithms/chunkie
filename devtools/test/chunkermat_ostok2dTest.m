% chunkermat_ostok2dTest0();

% function chunkermat_ostok2dTest0()

%CHUNKERMAT_OSTOK2DTEST
%
% test the matrix builder and do a basic solve

iseed = 8675309;
rng(iseed);

cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 1;
pref = []; 
pref.k = 16;
narms = 3;
amp = 0.25;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)
%%
% sources

ns = 10;
ts = 0.0+2*pi*rand(ns,1);
sources = starfish(ts,narms,amp);
sources = (1.5 + rand(1,ns)).*sources;
strengths = randn(2*ns,1);
sources_n = rand(2,ns);

% targets

nt = 100;
ts = 0.0+2*pi*rand(nt,1);
targets = starfish(ts,narms,amp);
targets = targets.*repmat(rand(1,nt),2,1)*0.8;

plot(chnkr, 'r.'); hold on;
plot(targets(1,:), targets(2,:), 'kx');
plot(sources(1,:), sources(2,:), 'bo');
axis equal

%% 
zk = 1.3;
kerns = kernel('ostok', 's', zk);

% eval u on bdry


srcinfo = []; srcinfo.r = sources; 
kernmats = kerns.eval(srcinfo, chnkr);
ubdry = kernmats*strengths;

%%
% eval u at targets

targinfo = []; targinfo.r = targets;
targets_n = rand(2, nt); 
targets_n = targets_n./sqrt(targets_n(1,:).^2+targets_n(2,:).^2);
targinfo.n = targets_n;
kernmatstarg = kerns.eval(srcinfo, targinfo);
utarg = kernmatstarg*strengths;
%%
% solve

% build Oscillatory Stokes dirichlet matrix
fkern = kernel('ostok', 's', zk);
start = tic; sys = chunkermat(chnkr, fkern);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

rhs = ubdry; rhs = rhs(:);
start = tic; sol = gmres(sys,rhs,[],1e-14,100); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

% evaluate at targets and compare

opts.usesmooth=false;
opts.verb=false;
start=tic; Dsol = chunkerkerneval(chnkr, fkern, sol, targets, opts); 
t1 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n',t1)

relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
fprintf('relative frobenius error %5.2e\n',relerr);

assert(relerr < 1e-10);

% end