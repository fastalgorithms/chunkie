
%CHUNKERMAT_HELM2DFMMTEST
%
% test the matrix builder and the fmm for post processing and do a basic solve

clearvars; close all;
iseed = 8675309;
rng(iseed);

addpaths_loc();

cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 1;
pref = []; 
pref.k = 32;
narms = 3;
amp = 0.25;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

zk = 1.1;

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
targets = targets.*repmat(rand(1,nt),2,1)*0.5;

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

kerns = @(s,t) chnk.helm2d.kern(zk,s,t,'s');

% eval u on bdry

targs = chnkr.r; targs = reshape(targs,2,chnkr.k*chnkr.nch);
targstau = tangents(chnkr); 
targstau = reshape(targstau,2,chnkr.k*chnkr.nch);

srcinfo = []; srcinfo.r = sources; 
targinfo = []; targinfo.r = targs;
kernmats = kerns(srcinfo,targinfo);
ubdry = kernmats*strengths;

% eval u at targets

targinfo = []; targinfo.r = targets;
kernmatstarg = kerns(srcinfo,targinfo);
utarg = kernmatstarg*strengths;

%

% build helmholtz dirichlet matrix

fkern = @(s,t) chnk.helm2d.kern(zk,s,t,'D');
start = tic; D = chunkermat(chnkr,fkern);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = -0.5*eye(chnkr.k*chnkr.nch) + D;

rhs = ubdry; rhs = rhs(:);
start = tic; sol = gmres(sys,rhs,[],1e-14,100); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

start = tic; sol2 = sys\rhs; t1 = toc(start);

fprintf('%5.2e s : time for dense backslash solve\n',t1)

err = norm(sol-sol2,'fro')/norm(sol2,'fro');

fprintf('difference between direct and iterative %5.2e\n',err)

% evaluate at targets and compare
wchnkr = chnkr.wts;

% evaluate at targets using FMM and compare
iffmm = 0;

if(iffmm)
    sol_use = sol2.*wchnkr(:);
    eps = 1e-6;
    pgt = 1;
    pot = chnk.helm2d.fmm(eps,zk,chnkr,targets,'D',sol_use);

    relerr = norm(utarg-pot,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
    relerr2 = norm(utarg-pot,'inf')/dot(abs(sol(:)),wchnkr(:));
    fprintf('relative frobenius error %5.2e\n',relerr);
    fprintf('relative l_inf/l_1 error %5.2e\n',relerr2);

    assert(relerr < eps);
end

