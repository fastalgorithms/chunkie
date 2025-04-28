chunkermat_stok_tractiontest0();


function chunkermat_stok_tractiontest0()

%CHUNKERMAT_HELM2DTEST
%
% test the matrix builder and do a basic solve

iseed = 8675309;
rng(iseed);


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
sources_n = rand(2,ns);
strengths = randn(2*ns,1);

% targets

nt = 10000;
ts = 0.0+2*pi*rand(nt,1);
targets = starfish(ts,narms,amp);
targets = targets.*repmat(rand(1,nt),2,1)*0.75;

plot(chnkr, 'r.'); hold on;
plot(targets(1,:), targets(2,:), 'kx')
axis equal
%
mu = 1.1; 
kerndt = kernel('stok', 'dtrac', mu);
kernd = kernel('stok', 'dvel', mu);
kerns = kernel('stok', 'svel', mu);

% eval u on bdry

targs = chnkr.r; targs = reshape(targs,2,chnkr.k*chnkr.nch);
srcinfo = []; srcinfo.r = sources; srcinfo.n = sources_n;

targinfo = []; targinfo.r = targs; targinfo.n = chnkr.n(:,:);
kernmats = kerndt.eval(srcinfo,targinfo);
ubdry = kernmats*strengths;

% eval u at targets

targinfo = []; targinfo.r = targets;
targets_n = rand(2, nt);
kernmatstarg = kernd.eval(srcinfo,targinfo);
utarg = kernmatstarg*strengths;


% build Stokes traction matrix
coefs = [1,1];
fkern = kernel('stok', 'strac', mu, coefs);
start = tic; D = chunkermat(chnkr, fkern);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = 0.5*eye(size(D,1)) + D;


rhs = ubdry; rhs = rhs(:);
start = tic; sol = gmres(sys,rhs,[],1e-14,100); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

% evaluate at targets and compare

opts.usesmooth=false;
opts.verb=false;
opts.forcefmm=true;
start=tic; Dsol = chunkerkerneval(chnkr,kerns,sol,targets,opts); 
t1 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n',t1)

tperp = zeros(2,nt);
tperp(1,:) = targets(2,:);
tperp(2,:) = -targets(1,:);

tperp = tperp(:);

A = [repmat(eye(2), [nt,1]), tperp];
rbd = A \ (utarg(:) - Dsol(:));

res = norm(A*rbd - (utarg(:) - Dsol(:)));

Dsol2 = reshape(Dsol, [2,nt]);
Dsol2(1,:) = Dsol2(1,:) + rbd(1) + rbd(3)*targets(2,:);
Dsol2(2,:) = Dsol2(2,:) + rbd(2) - rbd(3)*targets(1,:);


relerr = norm(utarg(:)-Dsol2(:),'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
fprintf('relative frobenius error %5.2e\n',relerr);

assert(relerr < 1e-10);


end


