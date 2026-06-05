% chunkermat_ostok2dTest0();

% function chunkermat_ostok2dTest0()

% CHUNKERMAT_OSTOK2DTEST
%
% test the matrix builder and do a basic solve
clear;

iseed = 8675309;
rng(iseed);

cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 1;
pref = []; 
pref.k = 20;
narms = 3;
amp = 0.25;
start = tic; 
chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1);

% sources

ns = 10;
ts = 0.0+2*pi*rand(ns,1);
sources = starfish(ts,narms,amp);
sources = 3.0*sources;
strengths = randn(2*ns,1);
sources_n = rand(2,ns);

% targets

nt = 100;
ts = 0.0+2*pi*rand(nt,1);
targets = starfish(ts,narms,amp);
targets = targets.*repmat(rand(1,nt),2,1)*0.8;

plot(chnkr, 'r.'); hold on;
plot(targets(1,:), targets(2,:), 'kx')
hold on;
plot(sources(1,:), sources(2,:), 'bo')
axis equal




targs = chnkr.r; targs = reshape(targs,2,chnkr.k*chnkr.nch);
targstau = tangents(chnkr); 
targstau = reshape(targstau,2,chnkr.k*chnkr.nch);

plot(chnkr, 'r.'); hold on;
plot(targets(1,:), targets(2,:), 'kx');
plot(sources(1,:), sources(2,:), 'bo');hold on;

t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)


zk = 0.3;
kerns = kernel('ostok', 's', zk);
kernd = kernel('ostok', 'd', zk);

% eval u on bdry


srcinfo = []; 
srcinfo.r = sources; 
srcinfo.n = sources_n;
kernmats = kerns.eval(srcinfo, chnkr); 
ubdry = kernmats*strengths; 

% eval u at targets

targinfo = []; 
targinfo.r = targets;
targets_n = rand(2, nt); 
targets_n = targets_n./sqrt(targets_n(1,:).^2+targets_n(2,:).^2);
targinfo.n = targets_n;
kernmatstargs = kerns.eval(srcinfo, targinfo); 

utarg = kernmatstargs*strengths; 


% solve

fkernd = kernel('ostok', 'd', zk);  
fkerns = kernel('ostok', 's', zk);
coefs = rand(1,2);
fkernc = kernel('ostok', 'c', zk, coefs);


start = tic; 
S = chunkermat(chnkr, fkerns);
D = chunkermat(chnkr, fkernd);
C = chunkermat(chnkr, fkernc);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sysd = -0.5*eye(size(D,1)) + D;
sysd = sysd + normonesmat(chnkr)/sum(chnkr.wts(:)); 

syss = S;

sysc = -0.5*coefs(1)*eye(size(D,1)) + C;
sysc = sysc + coefs(1)*normonesmat(chnkr)/sum(chnkr.wts(:));


rhs = ubdry; 
rhs = rhs(:);

start = tic; 
sold = gmres(sysd,rhs,[],1e-12,1000); 
sols = gmres(syss,rhs,[],1e-12,1000); 
solc = gmres(sysc,rhs,[],1e-12,1000); 


t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1);

% evaluate at targets and compare

opts.usesmooth=false;
opts.verb=false;

% SINGLE LAYER TEST 

fkerns =  kernel('ostok', 's', zk);
Ssol = chunkerkerneval(chnkr, fkerns, sols, targets, opts);
relerr = norm(utarg-Ssol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
fprintf('S_ relative frobenius error %5.2e\n', relerr);


% DOUBLE LAYER TEST

fkernd =  kernel('ostok', 'd', zk);
Dsol = chunkerkerneval(chnkr, fkernd, sold, targets, opts);
relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
fprintf('D_relative frobenius error %5.2e\n', relerr);

% COMBINED FIELD TEST

fkernc =  kernel('ostok', 'c', zk, coefs);
Dsol = chunkerkerneval(chnkr, fkernc, solc, targets, opts);
relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
fprintf('C_relative frobenius error %5.2e\n', relerr);