
%FLAMPROXYBYLEVELTEST
%
% test the FLAM matrix builder with level-dependent proxy points

clearvars; close all;
iseed = 8675309;
rng(iseed);

zk  = 50;
tol = 1e-14;
proxybylevel = true;

% addpaths_loc();

cparams = [];
cparams.maxchunklen = 4.0/abs(zk);
cparams.nover = 0;
chnkr = chunkerfunc(@(t) starfish(t), cparams, []);
chnkr = sort(chnkr);
s = chnkr.arclengthfun;

npt = chnkr.npt;
compute_dense = (npt < 5000);

wts = weights(chnkr); wts = wts(:);

% sources

ns = 100;
ts = 0.0+2*pi*rand(ns,1);
rs = (3.0 + 3*rand(ns,1));
sources = rs' .* [cos(ts)';sin(ts)'];
strengths = randn(ns,1);

% targets

nt = 1;
ts = 0.0+2*pi*rand(nt,1);
targets = 0.2*[cos(ts)'; sin(ts)'];
targets = targets.*repmat(rand(1,nt),2,1);

% plot geo and sources

figure(1)
clf
hold off
plot(chnkr)
hold on
scatter(sources(1,:),sources(2,:),'o')
scatter(targets(1,:),targets(2,:),'x')
axis equal 

% declare kernels

kerns = kernel('helm','s', zk);
kernd = -2*kernel('helm','d',zk);

% eval u on bdry

srcinfo = []; srcinfo.r = sources;
targinfo = []; targinfo.r = chnkr.r(:,:); 
targinfo.d = chnkr.d(:,:);
ubdry = kerns.fmm(tol,srcinfo,targinfo,strengths);

% eval u at targets

targinfo = []; targinfo.r = targets;
utarg = kerns.fmm(tol,srcinfo,targinfo,strengths);

%%

% build laplace dirichlet matrix * (-2)
if compute_dense
    start = tic; 
    D = chunkermat(chnkr,kernd);
    t1 = toc(start);
    
    fprintf('%5.2e s : time to assemble matrix\n',t1)
    
    sys = eye(npt) + D;
end

opts = [];
opts.rank_or_tol = tol;
opts.proxybylevel = proxybylevel;
opts.flamtype = 'rskelf';
opts.verb = true;

start = tic;
F = chunkerflam(chnkr,kernd,1.0,opts);
t1 = toc(start);

fprintf('%5.2e s : time for FLAM rskelf compress\n',t1)

rhs = ubdry;
start = tic; sig_rskelf = rskelf_sv(F,rhs); t1 = toc(start);

fprintf('%5.2e s : time for rskelf solve \n',t1)

%%

opts = [];
opts.rank_or_tol = tol;
opts.proxybylevel = proxybylevel;
opts.flamtype = 'rskel';
opts.verb = true;

start = tic;
F2 = chunkerflam(chnkr,kernd,1.0,opts);
t1 = toc(start);

fprintf('%5.2e s : time for FLAM rskel compress\n',t1)

rhs = ubdry;
start = tic; 
mtvc_rskel = rskel_mv(F2,rhs); 
t1 = toc(start);

fprintf('%5.2e s : time for rskel apply \n',t1)

if compute_dense
    start = tic;
    mtvc_dense = sys * rhs;
    t1 = toc(start);

    fprintf('%5.2e s : time for dense apply\n',t1)

    err = norm(mtvc_dense-mtvc_rskel,'fro')/norm(mtvc_dense,'fro');
        
    fprintf('difference between rskel and dense applies %5.2e\n\n',err)
    
    assert(err < 10*tol);

    start = tic; 
    sig_dense = gmres(sys,rhs,[],1e-14,min(npt, 1000)); 
    t1 = toc(start);

    fprintf('%5.2e s : time for dense gmres\n',t1)

    err = norm(sig_dense-sig_rskelf,'fro')/norm(sig_dense,'fro');
    
    fprintf('difference between rskelf and dense GMRES densities %5.2e\n\n',err)
    
    assert(err < 1e-12 || err < 10*tol);
end

%%

% evaluate at targets using adaptive quadrature

opts = [];
opts.eps = tol;
opts.forceadap = true;
opts.verb = true;

start=tic; 
sol_adap = chunkerkerneval(chnkr,kernd,sig_rskelf,targets,opts); 
t1 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n',t1)

relerr = norm(utarg-sol_adap,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
relerr2 = norm(utarg-sol_adap,'inf')/dot(abs(sig_rskelf(:)),wts(:));
fprintf('relative frobenius error %5.2e\n',relerr);
fprintf('relative l_inf/l_1 error %5.2e\n',relerr2);

%% evaluate at targets using FLAM's ifmm

opts = [];
opts.eps = tol;
opts.accel = true;
opts.verb = true;
opts.proxybylevel = proxybylevel;

% remove fmm to force FLAM ifmm usage
kernifmm = kernd;
kernifmm.fmm = [];

start=tic; 
sol_ifmm = chunkerkerneval(chnkr,kernifmm,sig_rskelf,targets,opts); 
t1 = toc(start);
fprintf('%5.2e s : time to eval at targs (ifmm)\n',t1)

relerr = norm(utarg-sol_ifmm,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
relerr2 = norm(utarg-sol_ifmm,'inf')/dot(abs(sig_rskelf(:)),wts(:));
fprintf('relative frobenius error %5.2e\n',relerr);
fprintf('relative l_inf/l_1 error %5.2e\n',relerr2);
