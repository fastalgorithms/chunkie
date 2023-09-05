
%CHUNKERMAT_AXISSYMHELM2DTEST
%
% test the matrix builder and do a basic solve

clearvars; close all;
iseed = 8675309;
rng(iseed);

addpaths_loc();

iftorus = 0;

cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 1;

if ~iftorus
    cparams.ta = -pi/2;
    cparams.tb = pi/2;
    ctr = [0 0];
    cparams.ifclosed = false;
else
    cparams.ta = 0;
    cparams.tb = 2*pi;
    ctr = [3 0];
    cparams.ifclosed = true;
end


cparams.maxchunklen = 0.5;
pref = []; 
pref.k = 16;
narms = 0;
amp = 0.25;


start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp,ctr),cparams,pref); 
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

zk = 80.1;

% sources

ns = 10;
ts = -pi/2+pi*rand(ns,1);
sources = starfish(ts,narms,amp);
sources = 0.5*sources;
sources(1,:) = sources(1,:) + ctr(1);
strengths = randn(ns,1);

% targets

nt = 100;
ts = -pi/2+pi*rand(nt,1);
targets = starfish(ts,narms,amp);
targets = targets .* (1 + 0.5*repmat(rand(1, nt), 2, 1));    
targets(1,:) = targets(1,:) + ctr(1);

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

kerns = @(s,t) chnk.axissymhelm2d.kern(zk,s,t, [0,0], 's');

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

% build laplace dirichlet matrix

fkern = kernel('axissymh', 'c', zk);
fkern = @(s,t) chnk.axissymhelm2d.kern(zk,s,t,[0 0],'c', [1 -1j]);
start = tic; D = chunkermat(chnkr,fkern);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = 0.5*eye(chnkr.k*chnkr.nch) + D;

rhs = ubdry; rhs = rhs(:);

start = tic; sol = gmres(sys,rhs,[],1e-14,200); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

start = tic; sol2 = sys\rhs; t1 = toc(start);

fprintf('%5.2e s : time for dense backslash solve\n',t1)

err = norm(sol-sol2,'fro')/norm(sol2,'fro');

fprintf('difference between direct and iterative %5.2e\n',err)

% evaluate at targets and compare

opts.usesmooth=false;
opts.verb=false;
opts.quadkgparams = {'RelTol',1e-16,'AbsTol',1.0e-16};
start=tic; Dsol = chunkerkerneval(chnkr,fkern,sol2,targets,opts); 
t1 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n',t1)

%

wchnkr = chnkr.wts;

relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
relerr2 = norm(utarg-Dsol,'inf')/dot(abs(sol(:)),wchnkr(:));
fprintf('relative frobenius error %5.2e\n',relerr);
fprintf('relative l_inf/l_1 error %5.2e\n',relerr2);

assert(relerr < 1e-7);
