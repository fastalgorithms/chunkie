
%FLAMUTILITIESTEST
%
% test the FLAM matrix builder and do a basic solve

clearvars; close all;
iseed = 8675309;
rng(iseed);

addpaths_loc();

cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 5;
pref = []; 
pref.k = 30;
narms = 3;
amp = 0.25;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

wts = weights(chnkr);

fprintf('%5.2e s : time to build geo\n',t1)

% sources

ns = 10;
ts = 0.0+2*pi*rand(ns,1);
sources = starfish(ts,narms,amp);
sources = 3.0*sources;
strengths = randn(ns,1);

% targets

nt = 1;
ts = 0.0+2*pi*rand(nt,1);
targets = starfish(ts,narms,amp);
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

%

kerns = @(s,t) chnk.lap2d.kern(s,t,'s');

% eval u on bdry

srcinfo = []; srcinfo.r = sources;
targinfo = []; targinfo.r = chnkr.r(:,:); 
targinfo.d = chnkr.d(:,:);
kernmats = kerns(srcinfo,targinfo);
ubdry = kernmats*strengths;

% eval u at targets

targinfo = []; targinfo.r = targets;
kernmatstarg = kerns(srcinfo,targinfo);
utarg = kernmatstarg*strengths;

%

% build laplace dirichlet matrix

fkern = @(s,t) chnk.lap2d.kern(s,t,'D');
start = tic; D = chunkermat(chnkr,fkern);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = -0.5*eye(chnkr.k*chnkr.nch) + D;

% build sparse tridiag part 

opts.nonsmoothonly = true;
start = tic; spmat = chunkermat(chnkr,fkern,opts);
t1 = toc(start);
fprintf('%5.2e s : time to build tridiag\n',t1)

spmat = spmat -0.5*speye(chnkr.k*chnkr.nch);

% test matrix entry evaluator
start = tic; 
opdims = [1 1];
sys2 = chnk.flam.kernbyindex(1:chnkr.npt,1:chnkr.npt,chnkr,wts,fkern,opdims,spmat);
t1 = toc(start);

fprintf('%5.2e s : time for mat entry eval on whole mat\n',t1)


xflam = chnkr.r(:,:);
matfun = @(i,j) chnk.flam.kernbyindex(i,j,chnkr,wts,fkern,opdims,spmat);
[pr,ptau,pw,pin] = chnk.flam.proxy_square_pts();
ifaddtrans = true;
pxyfun = @(x,slf,nbr,l,ctr) chnk.flam.proxyfun(slf,nbr,l,ctr,chnkr,wts, ...
    fkern,opdims,pr,ptau,pw,pin,ifaddtrans);


start = tic; F = rskelf(matfun,xflam,200,1e-14,pxyfun); t1 = toc(start);
F = chunkerflam(chnkr,fkern,-0.5);

fprintf('%5.2e s : time for flam rskelf compress\n',t1)

pxyfunr = @(rc,rx,cx,slf,nbr,l,ctr) chnk.flam.proxyfunr(rc,rx,slf,nbr,l, ...
        ctr,chnkr,wts,fkern,opdims,pr,ptau,pw,pin);

opts = [];
start = tic; F2 = rskel(matfun,xflam,xflam,200,1e-14,pxyfunr,opts); t1 = toc(start);

fprintf('%5.2e s : time for flam rskel compress\n',t1)

afun = @(x) rskelf_mv(F,x);



err2 = norm(sys2-sys,'fro')/norm(sys,'fro');
fprintf('%5.2e   : fro error of build \n',err2)

rhs = ubdry; rhs = rhs(:);
start = tic; sol = gmres(sys,rhs,[],1e-14,100); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

rhs = ubdry; rhs = rhs(:);
start = tic; sol3 = rskelf_sv(F,rhs); t1 = toc(start);

fprintf('%5.2e s : time for rskelf_sv \n',t1)

err = norm(sol-sol3,'fro')/norm(sol,'fro');

fprintf('difference between fast-direct and iterative %5.2e\n',err)

assert(err < 1e-10);

% evaluate at targets and compare

opts.usesmooth=false;
opts.verb=false;
opts.quadkgparams = {'RelTol',1e-16,'AbsTol',1.0e-16};
start=tic; Dsol = chunkerkerneval(chnkr,fkern,sol3,targets,opts); 
t1 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n',t1)

%

wchnkr = weights(chnkr);

relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
relerr2 = norm(utarg-Dsol,'inf')/dot(abs(sol(:)),wchnkr(:));
fprintf('relative frobenius error %5.2e\n',relerr);
fprintf('relative l_inf/l_1 error %5.2e\n',relerr2);

