chunkermat_flex2d_exteriorfreeTest0();


function chunkermat_flex2d_exteriorfreeTest0()

%CHUNKERMAT_FLEX2D_EXTERIORFREETEST
%
% test the matrix builder and do a basic solve

iseed = 8675309;
rng(iseed);

% PDE coefficients: (a \Delta^2 - b \Delta - c) u = 0
a = 1.1;
b = 0.7;
c = 1/pi;
nu = 0.3;

zk1 = sqrt((- b + sqrt(b^2 + 4*a*c)) / (2*a));
zk2 = sqrt((- b - sqrt(b^2 + 4*a*c)) / (2*a));

zk = [zk1 zk2];

% zk = 3;

cparams = [];
% cparams.eps = 1.0e-6;
cparams.nover = 1;
cparams.maxchunklen = 4.0/max(abs(zk));
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

kern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 's');
kern2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_bcs',nu);

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

% defining free plate kernels

fkern1 =  @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate', nu);        % build the desired kernel
double = @(s,t) chnk.lap2d.kern(s,t,'d');
hilbert = @(s,t) chnk.lap2d.kern(s,t,'hilb');

opts = [];
opts.sing = 'log';

opts2 = [];
opts2.sing = 'pv';

% building system matrix

start = tic;
sysmat1 = chunkermat(chnkr,fkern1, opts);
D = chunkermat(chnkr, double, opts);
H = chunkermat(chnkr, hilbert, opts2);     

sysmat = zeros(2*chnkr.npt);
sysmat(1:2:end,1:2:end) = sysmat1(1:4:end,1:2:end) + sysmat1(3:4:end,1:2:end)*H  - 2*((1+nu)/2)^2*D*D;
sysmat(2:2:end,1:2:end) = sysmat1(2:4:end,1:2:end) + sysmat1(4:4:end,1:2:end)*H;
sysmat(1:2:end,2:2:end) = sysmat1(1:4:end,2:2:end) + sysmat1(3:4:end,2:2:end);
sysmat(2:2:end,2:2:end) = sysmat1(2:4:end,2:2:end) + sysmat1(4:4:end,2:2:end);

D = [-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2];  % jump matrix 
D = kron(eye(chnkr.npt), D);

sys =  D + sysmat;
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

ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_eval', nu); 

dens_comb = zeros(3*chnkr.npt,1);
dens_comb(1:3:end) = sol(1:2:end);
dens_comb(2:3:end) = H*sol(1:2:end);
dens_comb(3:3:end) = sol(2:2:end);

start1 = tic;
Dsol = chunkerkerneval(chnkr, ikern,dens_comb,targets);
t2 = toc(start1);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n',t2)

% calculate error

wchnkr = chnkr.wts;

relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
relerr2 = norm(utarg-Dsol,'inf')/dot(abs(sol(1:2:end)+sol(2:2:end)),wchnkr(:));
fprintf('relative frobenius error %5.2e\n',relerr);
fprintf('relative l_inf/l_1 error %5.2e\n',relerr2);

assert(relerr < 1e-9);

end


