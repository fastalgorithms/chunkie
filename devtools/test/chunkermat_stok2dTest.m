
%CHUNKERMAT_HELM2DTEST
%
% test the matrix builder and do a basic solve

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
strengths = randn(2*ns,1);
sources_n = rand(2,ns);

% targets

nt = 10000;
ts = 0.0+2*pi*rand(nt,1);
targets = starfish(ts,narms,amp);
targets = targets.*repmat(rand(1,nt),2,1)*0.75;
% targets = targets.*repmat(rand(1,nt),2,1)*0.99;

plot(chnkr, 'r.'); hold on;
plot(targets(1,:), targets(2,:), 'kx')
axis equal
%
mu = 1.3; 
kerns = kernel('stok', 'd', mu);

% eval u on bdry

targs = chnkr.r; targs = reshape(targs,2,chnkr.k*chnkr.nch);
targstau = tangents(chnkr); 
targstau = reshape(targstau,2,chnkr.k*chnkr.nch);

srcinfo = []; srcinfo.r = sources; 
srcinfo.n = sources_n;
targinfo = []; targinfo.r = targs;
kernmats = kerns.eval(srcinfo,targinfo);
ubdry = kernmats*strengths;

% eval u at targets

targinfo = []; targinfo.r = targets;
targets_n = rand(2, nt);
kernmatstarg = kerns.eval(srcinfo,targinfo);
utarg = kernmatstarg*strengths;

% solve

% build Stokes dirichlet matrix
coefs = [1,1];
fkern = kernel('stok', 'cvel', mu, coefs);
start = tic; D = chunkermat(chnkr,fkern);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = -0.5*eye(size(D,1)) + D;

sys = sys + normonesmat(chnkr)/sum(chnkr.wts(:));

rhs = ubdry; rhs = rhs(:);
start = tic; sol = gmres(sys,rhs,[],1e-14,100); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

% evaluate at targets and compare

opts.usesmooth=false;
opts.verb=false;
opts.forcefmm=true;
start=tic; Dsol = chunkerkerneval(chnkr,fkern,sol,targets,opts); 
t1 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n',t1)

% just test stokes slp velocity pquad...
% targets = targets.*repmat(rand(1,nt),2,1)*0.99 seems to break below stokes kernel test...
fkerns = kernel('stok', 'svel', mu);
Ssol = chunkerkerneval(chnkr,fkerns,sol,targets,opts); 
opts.forcepquad=true;
opts.side = 'i';
Ssol_pquad = chunkerkerneval(chnkr,fkerns,sol,targets,opts); 
err = abs(Ssol - Ssol_pquad);
assert(norm(err) < 1e-10);
% figure(1),clf,
% plot(chnkr), hold on
% scatter(targets(1,:),targets(2,:),[],log10(err(1:2:end))); colorbar
opts.forcepquad=false;
%

wchnkr = chnkr.wts;

relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
fprintf('relative frobenius error %5.2e\n',relerr);

assert(relerr < 1e-10);

% Test pressure in the bulk
tol = 1e-11;
ww = repmat(wchnkr(:).', [2,1]);
ww = ww(:);
soluse = sol.*ww;
[~, p] = fkern.fmm(tol, chnkr, targets, soluse);
fkernp = kernel('stok', 'dpres', mu, coefs);

pex = fkernp.eval(srcinfo, targinfo)*strengths;


p = p - p(1);
pex = pex - pex(1);
relerr = norm(p - pex)/norm(pex);
assert(relerr < 1e-10)

% Test pressure in a different way now
fkernp_c = kernel('stokes', 'cpres', mu, coefs);
p = chunkerkerneval(chnkr,fkernp_c,sol,targets,opts); 
p = p - p(1);
relerr = norm(p - pex)/norm(pex);

assert(relerr < 1e-10)


% Test gradient in the bulk
tol = 1e-11;
ww = repmat(wchnkr(:).', [2,1]);
ww = ww(:);
soluse = sol.*ww;
[~, ~, g] = fkern.fmm(tol, chnkr, targets, soluse);

fkerng = kernel('stok', 'dgrad', mu, coefs);
gex = fkerng.eval(srcinfo, targinfo)*strengths;

g = reshape(g, [2,2,nt]);
gex = reshape(gex, [2,2,nt]);

relerr = norm(g(:) - gex(:))/norm(gex(:));
assert(relerr < 1e-10)

fkerng_c = kernel('stokes', 'cgrad', mu, coefs);
g = chunkerkerneval(chnkr,fkerng_c, sol, targets, opts); 
g = reshape(g, [2,2,nt]);
relerr = norm(g(:) - gex(:))/norm(gex(:));

assert(relerr < 1e-10)

