chunkermat_biharmonic_general_supportedTest0();


function chunkermat_biharmonic_general_supportedTest0()

%CHUNKERMAT_BIHARMONIC_GENERAL_SUPPORTEDTEST
%
% test the matrix builder and do a basic solve

iseed = 8675309;
rng(iseed);

a = 1.7;
b = 0.8;
c = 1/pi;
nu = 0.3;

zk1 = sqrt((- b + sqrt(b^2 + 4*a*c)) / (2*a));
zk2 = sqrt((- b - sqrt(b^2 + 4*a*c)) / (2*a));

cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 1;
cparams.maxchunklen = 4.0/zk1;
pref = []; 
pref.k = 16;
narms = 3;
amp = 0.25;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
chnkr = makedatarows(chnkr,2);
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
kern2 = @(s,t) chnk.flex2d.kern([], s, t, 'supported_plate_bcs_general',a,b,c,nu);

% eval boundary conditions on bdry

srcinfo = []; srcinfo.r = sources; 
targinfo = chnkr;

ubdry = kern2(srcinfo,targinfo);
rhs = ubdry*strengths;

% eval u at targets

srcinfo = []; srcinfo.r = sources; 
targinfo = []; targinfo.r = targets;
kernmatstarg = kern1(srcinfo,targinfo);
utarg = kernmatstarg*strengths;

% calculating curvature info

kappa = signed_curvature(chnkr);
kp = arclengthder(chnkr,kappa);
kpp = arclengthder(chnkr,kp);

chnkr.data(1,:,:) = kp;
chnkr.data(2,:,:) = kpp;

% defining supported plate kernels

fkern1 =  @(s,t) chnk.flex2d.kern([], s, t, 'supported_plate_general_log',a,b,c,nu);           % build the desired kernel
fkern2 =  @(s,t) chnk.flex2d.kern([], s, t, 'supported_plate_smooth',nu);           % build the desired kernel

opts1 = [];
opts1.sing = 'log';

opts2 = [];
opts2.quad = 'native';
opts2.sing = 'smooth';

start = tic;
M = chunkermat(chnkr,fkern1, opts1);
M2 = chunkermat(chnkr,fkern2, opts2);

c0 = (nu - 1)*(nu + 3)*(2*nu - 1)/(2*(3 - nu));

M(2:2:end,1:2:end) = M(2:2:end,1:2:end) + M2 + c0.*kappa(:).^2.*eye(chnkr.npt) - b/(2*a)*eye(chnkr.npt); % extra term shows up for the general problem
M = M - 0.5*eye(2*chnkr.npt);

sys =  M;
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

% solve linear system

start = tic; sol = gmres(sys,rhs,[],1e-10,200); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

start = tic; sol2 = sys\rhs; t1 = toc(start);

fprintf('%5.2e s : time for dense backslash solve\n',t1)

err = norm(sol-sol2,'fro')/norm(sol2,'fro');

fprintf('difference between direct and iterative %5.2e\n',err)

% evaluate at targets and compare

ikern = @(s,t) chnk.flex2d.kern([], s, t, 'supported_plate_general_eval',a,b,c,nu);                              % build the kernel of evaluation          

start1 = tic;
Dsol = chunkerkerneval(chnkr, ikern, sol, targets);t2 = toc(start1);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n',t2)

% calculate error

wchnkr = chnkr.wts;

relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
relerr2 = norm(utarg-Dsol,'inf')/dot(abs(sol(1:2:end)+sol(2:2:end)),wchnkr(:));
fprintf('relative frobenius error %5.2e\n',relerr);
fprintf('relative l_inf/l_1 error %5.2e\n',relerr2);

assert(relerr < 1e-9);

end


