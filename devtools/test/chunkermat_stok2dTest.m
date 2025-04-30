chunkermat_stok2dTest0();


function chunkermat_stok2dTest0()

%CHUNKERMAT_HELM2DTEST
%
% test the matrix builder and do a basic solve

iseed = 8675309;
rng(iseed);


cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 1;
pref = []; 
pref.k = 20;
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
targets = targets.*repmat(rand(1,nt),2,1)*0.8;
% targets = targets.*repmat(rand(1,nt),2,1)*0.99;

plot(chnkr, 'r.'); hold on;
plot(targets(1,:), targets(2,:), 'kx')
axis equal
%
mu = 1.3;
% mu = 1;
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
targets_n = targets_n./sqrt(targets_n(1,:).^2+targets_n(2,:).^2);
targinfo.n = targets_n;
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

% test Stokes SLP velocity
fkerns = kernel('stok', 'svel', mu);
Ssol = chunkerkerneval(chnkr,fkerns,sol,targets,opts); 
Ssys = chunkerkernevalmat(chnkr,fkerns,targets,opts); 
err = abs(Ssol - Ssys*sol);
assert(norm(err) < 1e-10);
%
opts.forcepquad=true;
opts.side = 'i';
Ssol_pquad = chunkerkerneval(chnkr,fkerns,sol,targets,opts); 
err = abs(Ssol - Ssol_pquad);
assert(norm(err) < 1e-10);
opts.forcepquad=false;
%

% test Stokes DLP velocity
fkernd = kernel('stok', 'dvel', mu);
Dsol = chunkerkerneval(chnkr,fkernd,sol,targets,opts); 
%
opts.forcepquad=true;
opts.side = 'i';
Dsol_pquad = chunkerkerneval(chnkr,fkernd,sol,targets,opts); 
err = abs(Dsol - Dsol_pquad);
assert(norm(err) < 1e-10);
opts.forcepquad=false;

% test Stokes SLP traction 
fkernstrac = kernel('stok', 'strac', mu);
Stracsol = chunkerkerneval(chnkr,fkernstrac,sol,targinfo,opts);
% 
opts.forcepquad=true;
opts.side = 'i';
Stracsol_pquad = chunkerkerneval(chnkr,fkernstrac,sol,targinfo,opts); 
err = abs(Stracsol - Stracsol_pquad);
assert(norm(err) < 1e-10);
opts.forcepquad=false;

% test Stokes DLP traction 
fkerndtrac = kernel('stok', 'dtrac', mu);
Dtracsol = chunkerkerneval(chnkr,fkerndtrac,sol,targinfo,opts);
%
opts.forcepquad=true;
opts.side = 'i';
Dtracsol_pquad = chunkerkerneval(chnkr,fkerndtrac,sol,targinfo,opts); 
err = abs(Dtracsol - Dtracsol_pquad);
assert(norm(err) < 1e-8);
opts.forcepquad=false;

% test Stokes SLP pressure
fkernspres = kernel('stok', 'spres', mu);
Spressol = chunkerkerneval(chnkr,fkernspres,sol,targinfo,opts);
%
opts.forcepquad=true;
opts.side = 'i';
Spressol_pquad = chunkerkerneval(chnkr,fkernspres,sol,targinfo,opts); 
err = abs(Spressol - Spressol_pquad);
assert(norm(err) < 1e-10);
opts.forcepquad=false;

% test Stokes DLP pressure
fkerndpres = kernel('stok', 'dpres', mu);
Dpressol = chunkerkerneval(chnkr,fkerndpres,sol,targinfo,opts);
%
opts.forcepquad=true;
opts.side = 'i';
Dpressol_pquad = chunkerkerneval(chnkr,fkerndpres,sol,targinfo,opts); 
err = abs(Dpressol - Dpressol_pquad);
assert(norm(err) < 1e-10);
opts.forcepquad=false;


wchnkr = chnkr.wts;

Dsol = chunkerkerneval(chnkr,fkern,sol,targets,opts); 
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
assert(relerr < 1.1e-10)

% Test pressure in a different way now
fkernp_c = kernel('stokes', 'cpres', mu, coefs);
p = chunkerkerneval(chnkr,fkernp_c,sol,targets,opts); 
p = p - p(1);
relerr = norm(p - pex)/norm(pex);

assert(relerr < 1.1e-10)


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



end


