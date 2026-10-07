swirlfishTest0();


function swirlfishTest0()
%SWIRLFISHTEST test the swirlfish curve: derivatives, reduction to the 
% starfish at zero twist, handedness, chunking, area, and Green's identity
% for Laplace on a domain that is not star-shaped

seed = 8675309;
rng(seed);

narms = 5;
amp = 0.5;
ctr = [0.7;-1.2];
phi = 0.4;
scale = 2;
twist = 1;
fcurve = @(t) swirlfish(t,narms,amp,twist,ctr,phi,scale);

t = 2*pi*rand(1,50);
h = 1e-5;
[~,d,d2] = fcurve(t);
[rp,dp] = fcurve(t+h);
[rm,dm] = fcurve(t-h);
errd = norm(d-(rp-rm)/(2*h),'fro')/norm(d,'fro');
errd2 = norm(d2-(dp-dm)/(2*h),'fro')/norm(d2,'fro');
fprintf('swirlfish finite difference errors: d %5.2e, d2 %5.2e\n', ...
    errd,errd2);
assert(errd < 1e-7);
assert(errd2 < 1e-7);

% zero twist reduces to the starfish.

[r0,d0,d20] = swirlfish(t,narms,amp,0,ctr,phi,scale);
[r1,d1,d21] = starfish(t,narms,amp,ctr,phi,scale);
assert(norm(r0-r1,'fro') < 1e-12);
assert(norm(d0-d1,'fro') < 1e-12);
assert(norm(d20-d21,'fro') < 1e-11);

% handedness: flipping the sign of twist reflects the curve across the 
% x-axis (with ctr = 0 and phi = 0)

[rpos,dpos,d2pos] = swirlfish(-t,narms,amp,twist);
[rneg,dneg,d2neg] = swirlfish(t,narms,amp,-twist);
flip = [1 0; 0 -1];
assert(norm(flip*rpos-rneg,'fro') < 1e-12);
assert(norm(-flip*dpos-dneg,'fro') < 1e-12);
assert(norm(flip*d2pos-d2neg,'fro') < 1e-11);

% Test that geometry information is resolved

cparams = [];
cparams.eps = 1.0e-14;
pref = []; 
pref.k = 16;
aex = pi*(1+amp^2/2)*scale^2;

chnkr = chunkerfunc(@(t) swirlfish(t,narms,amp,twist,ctr,phi,scale), ...
        cparams,pref); 

kappa = signed_curvature(chnkr); kappa = kappa(:);

[~,~,v2c] = lege.exps(16);
kcoefs = v2c*reshape(kappa,[16,chnkr.nch]);
kcoefs = kcoefs ./ kcoefs(1,:);
kcoefs = abs(kcoefs);


figure(1); clf 
plot(chnkr); hold on;
quiver(chnkr)

figure(2); clf 
plot(log10(kcoefs))

assert(max(kcoefs(end,:)) < 1e-9)

% Test the derivatives agree

D = lege.dermat(16);        

r = chnkr.r;  d = chnkr.d;  d2 = chnkr.d2;
du  = zeros(size(r));  d2u = zeros(size(r));
for i = 1:chnkr.nch
    du(:,:,i)  = (D*r(:,:,i).').';
    d2u(:,:,i) = (D*d(:,:,i).').';
end
err_d  = norm(du(:)-d(:))/norm(d(:));
err_d2 = norm(d2u(:)-d2(:))/norm(d2(:));

drds = [reshape(arclengthder(chnkr,r(1,:,:)),1,[]);
        reshape(arclengthder(chnkr,r(2,:,:)),1,[])];
tau  = tangents(chnkr);
err_tau = norm(drds - tau(:,:),'fro')/sqrt(chnkr.npt);

Ds = diffmat(chnkr);                      
rss = (Ds*(Ds*r(:,:).')).';
kap = signed_curvature(chnkr);
nu  = [-tau(2,:); tau(1,:)];
err_kap = norm(rss - kap(:).'.*nu,'fro')/norm(rss,'fro');

assert(err_d < 1e-10)
assert(err_d2 < 1e-10)
assert(err_tau < 1e-10)
assert(err_kap < 1e-8)

% chunk up and check area, which the twist map preserves

narms = 3; cparams.eps = 1e-12;
for twist = [0.5 1.0 3.0]
    start = tic; 
    chnkr = chunkerfunc(@(t) swirlfish(t,narms,amp,twist,ctr,phi,scale), ...
        cparams,pref); 
    t1 = toc(start);

    [~,~,info] = sortinfo(chnkr);
    assert(info.ier == 0,'adjacency issues after chunk build swirlfish');

    err = abs(area(chnkr)-aex)/aex;
    fprintf(['%5.2e seconds to chunk swirlfish (twist %3.1f) with %d ' ...
        'chunks, area error %5.2e\n'],t1,twist,chnkr.nch,err);
    assert(err < 1e-10,'area of swirlfish should match starfish area');
end

% Green's identity for Laplace on the last (most curled) domain

ns = 50;
ts = 2*pi*rand(1,ns);
sources = ctr + 3*(1+amp)*scale*[cos(ts); sin(ts)];
strengths = randn(ns,1);

rmin = min(chnkr); rmax = max(chnkr);
targs = rmin + (rmax-rmin).*rand(2,2000);
in = chunkerinterior(chnkr,targs);
targets = targs(:,in);
assert(nnz(in) > 100,'expected more interior targets');

kernd = kernel('lap','d');
kerns = kernel('lap','s');
kernsprime = kernel('lap','sprime');

srcinfo = []; srcinfo.r = sources; 
targinfo = []; targinfo.r = chnkr.r(:,:); targinfo.d = chnkr.d(:,:);
densu = kerns.eval(srcinfo,targinfo)*strengths;
densun = kernsprime.eval(srcinfo,targinfo)*strengths;

targinfo = []; targinfo.r = targets;
utarg = kerns.eval(srcinfo,targinfo)*strengths;

Du = chunkerkerneval(chnkr,kernd,densu,targets);
Sun = chunkerkerneval(chnkr,kerns,densun,targets);
relerr = norm(utarg-(Sun-Du))/norm(utarg);
fprintf('swirlfish Green''s identity relative error %5.2e\n',relerr);
assert(relerr < 1e-10,'Green''s identity failed on swirlfish');

% Solve exterior free plate problem on a swirlfish

% PDE coefficients: (a \Delta^2 - b \Delta - c) u = 0
a = 1.1;
b = 0.7;
c = 1/pi;
nu = 0.3;

zk1 = sqrt((- b + sqrt(b^2 + 4*a*c)) / (2*a));
zk2 = sqrt((- b - sqrt(b^2 + 4*a*c)) / (2*a));

zk = [zk1 zk2];

% zk = 3;

narms = 3;
amp = 0.35;
ctr = [0;0];
phi = 0;
scale = 1.7;
twist = 1.2;

cparams = [];
% cparams.eps = 1.0e-6;
cparams.nover = 1;
cparams.maxchunklen = 4.0/max(abs(zk));
pref = []; 
pref.k = 16;
start = tic; chnkr = chunkerfunc(@(t) swirlfish(t,narms,amp,twist,ctr,phi,scale),cparams,pref); 
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

% targets

nt = 10;
ts = 0.0+2*pi*rand(nt,1);
targets = starfish(ts,narms,amp) + 0.5*[cos(ts.');sin(ts.')];
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

start = tic; sol = sys \ rhs; t1 = toc(start);

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

assert(relerr < 1e-8);


end
