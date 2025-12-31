chunkermat_biharmonic_general_clampedTest0();

function chunkermat_biharmonic_general_clampedTest0()

%CHUNKERMAT_BIHARMONIC_GENERALCLAMPEDTEST
%
% test the matrix builder and do a basic solve

iseed = 8675309;
rng(iseed);

a = 0.77;
b = 1.3;
c = 1/pi;

zk1 = sqrt((- b + sqrt(b^2 + 4*a*c)) / (2*a));
zk2 = sqrt((- b - sqrt(b^2 + 4*a*c)) / (2*a));

cparams = [];
% cparams.eps = 1.0e-12;
cparams.nover = 1;
cparams.maxchunklen = 4.0/zk1;
pref = []; 
pref.k = 16;
narms = 3;
amp = 0.25;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

% targets

nt = 10;
ts = 0.0+2*pi*rand(nt,1);
targets = starfish(ts,narms,amp);
targets = 3.0*targets;

% sources

ns = 3;
ts = 0.0+2*pi*rand(ns,1);
sources = starfish(ts,narms,amp);
sources = sources.*repmat(rand(1,ns),2,1);
strengths = randn(ns,1);

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

% defining kernels for rhs and analytic sol test

kern1 = @(s,t) chnk.flex2d.kern([], s, t, 's_general',a,b,c);
kern2 = @(s,t) chnk.flex2d.kern([], s, t, 'clamped_plate_bcs_general',a,b,c);

% eval u and du/dn on bdry

srcinfo = []; srcinfo.r = sources; 
targinfo = chnkr;

ubdry = kern2(srcinfo,targinfo);
rhs = ubdry*strengths;

% eval u at targets

targinfo = []; targinfo.r = targets;
kernmatstarg = kern1(srcinfo,targinfo);
utarg = kernmatstarg*strengths;

% build clamped plate matrix

fkern =  @(s,t) chnk.flex2d.kern([], s, t, 'clamped_plate_general',a,b,c);

kappa = signed_curvature(chnkr);
kappa = kappa(:);

opts = [];
opts.sing = 'log';

start = tic;
sys = chunkermat(chnkr,fkern, opts);
sys = sys - 0.5*eye(2*chnkr.npt);
sys(2:2:end,1:2:end) = sys(2:2:end,1:2:end) + kappa.*eye(chnkr.npt);

t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

start = tic; sol = gmres(sys,rhs,[],1e-10,200); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

start = tic; sol2 = sys\rhs; t1 = toc(start);

fprintf('%5.2e s : time for dense backslash solve\n',t1)

err = norm(sol-sol2,'fro')/norm(sol2,'fro');

fprintf('difference between direct and iterative %5.2e\n',err)

% evaluate at targets and compare

ikern = @(s,t) chnk.flex2d.kern([], s, t, 'clamped_plate_general_eval',a,b,c); 

start=tic; Dsol = chunkerkerneval(chnkr, ikern,sol,targets);
t1 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n',t1)

%

wchnkr = chnkr.wts;

relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
relerr2 = norm(utarg-Dsol,'inf')/dot(abs(sol(1:2:end)+sol(2:2:end)),wchnkr(:));
fprintf('relative frobenius error %5.2e\n',relerr);
fprintf('relative l_inf/l_1 error %5.2e\n',relerr2);

assert(relerr < 1e-9);
assert(relerr2 < 1e-9);


end


