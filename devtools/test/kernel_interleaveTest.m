clearvars; close all;
iseed = 8675309;
rng(iseed);

addpaths_loc();

% dumb test 

K0 = [kernel.lap2d('d'), kernel.nans()];
didfail = false;
try 
    K0 = interleave(K0);
catch
    didfail = true;
end
assert(didfail)

% solver test

zk = 1.1;

cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 1;
pref = []; 
pref.k = 16;
narms = 3;
amp = 0.25;
start = tic; chnkr = chunkerfunc(@(t) starfish(t, narms, amp), cparams, pref); 
t1 = toc(start);
chnkr = sort(chnkr);
wts = chnkr.wts; wts = wts(:);

l2scale = true;

fprintf('%5.2e s : time to build geo\n',t1)

% sources

ns = 10;
ts = 2*pi*rand(ns, 1);
sources = starfish(ts, narms, amp);
sources = 0.5*sources;
strengths = randn(ns, 1);

% targets

nt = 1000;
ts = 2*pi*rand(nt, 1);
targets = starfish(ts, narms, amp);
targets = targets .* (1 + 3*repmat(rand(1, nt), 2, 1));

% Plot everything

figure(1)
clf
hold off
plot(chnkr)
hold on
scatter(sources(1,:), sources(2,:), 'o')
scatter(targets(1,:), targets(2,:), 'x')
axis equal

% For solving exterior the Neumann boundary value problem, we use the
% following integral equation
%
% u = \beta(S_{k} + i \alpha D_{k} S_{ik}) \sigma
%
% with \beta = -1.0/(0.5 + 1i*0.25*alpha)
% 
% On imposing the boundary conditions, we get the following integral 
% equation
%
%  du/dn = I + \beta S_{k}'[\sigma] + 
%          i\beta \alpha(D'_{k} - D'_{ik})(S_{ik}[\sigma]) + 
%          i\beta \alpha(S_{ik}')^2 [\sigma];
% 
% Setting -S_{ik}[\sigma] + \tau = 0, and - S_{ik}'[\sigma] + \mu=0, 
% we have the following system of integral equations

% Set up kernels
alpha = 1;
c1 = -1/(0.5 + 1i*alpha*0.25);
c2 = -1i*alpha/(0.5 + 1i*alpha*0.25);
c3 = -1;
Sik    = kernel('helm', 's', 1i*zk);
Sikp   = kernel('helm', 'sprime', 1i*zk);
Skp    = kernel('helm', 'sprime', zk);
Sk     = kernel('helm', 's', zk);
Dk     = kernel('helm', 'd', zk);
Dkdiff = kernel('helmdiff', 'dprime', [zk 1i*zk]);
Z = kernel.zeros();
K = [ c1*Skp  c2*Dkdiff c2*Sikp ;
      c3*Sik  Z        Z        ;
      c3*Sikp Z        Z        ];
K = kernel(K);

% Set up boundary data

srcinfo  = []; srcinfo.r = sources; 
targinfo = []; targinfo.r = chnkr.r(:,:); targinfo.n = chnkr.n(:,:);
kernmats = Skp.eval(srcinfo, targinfo);
ubdry = kernmats*strengths;

npts = chnkr.npt;
nsys = 3*npts;
rhs = zeros(nsys, 1);


if(l2scale)
    rhs(1:3:end) = ubdry.*sqrt(wts);
else
    rhs(1:3:end) = ubdry;
end

% Form matrix
opts = [];
opts.l2scale = l2scale;
A = chunkermat(chnkr, K, opts) + eye(nsys);
start = tic;
sol = gmres(A, rhs, [], 1e-14, 100);
t1 = toc(start);

% Compute exact solution
srcinfo  = []; srcinfo.r  = sources;
targinfo = []; targinfo.r = targets;
kernmatstarg = Sk.eval(srcinfo, targinfo);
utarg = kernmatstarg*strengths;

% Compute solution using chunkerkerneval
% evaluate at targets and compare
Keval = c1*kernel([Sk 1i*alpha*Dk Z]);
opts.usesmooth = false;
opts.verb = false;
opts.quadkgparams = {'RelTol', 1e-16, 'AbsTol', 1.0e-16};

if(l2scale)
    wts_rep = repmat(wts(:).', K.opdims(1),1);
    wts_rep = wts_rep(:);
    sol = sol./sqrt(wts_rep);
end

start = tic;
Dsol = chunkerkerneval(chnkr, Keval, sol, targets, opts);
t2 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n', t2)


wchnkr = chnkr.wts;
wchnkr = repmat(wchnkr(:).', 3, 1);
relerr  = norm(utarg-Dsol) / (sqrt(chnkr.nch)*norm(utarg));
relerr2 = norm(utarg-Dsol, 'inf') / dot(abs(sol(:)), wchnkr(:));
fprintf('relative frobenius error %5.2e\n', relerr);
fprintf('relative l_inf/l_1 error %5.2e\n', relerr2);



% Test fast direct solver interfaces

% build sparse tridiag part 
opts.nonsmoothonly = true;
opts.rcip = true;
start = tic; spmat = chunkermat(chnkr, K, opts); t1 = toc(start);
fprintf('%5.2e s : time to build tridiag\n',t1)

spmat = spmat + speye(nsys);

% test matrix entry evaluator
start = tic; 
opdims = K.opdims;
sys2 = chnk.flam.kernbyindex(1:nsys, 1:nsys, chnkr, K, opdims, ...
    spmat, l2scale);


t1 = toc(start);

fprintf('%5.2e s : time for mat entry eval on whole mat\n',t1)

err2 = norm(sys2-A,'fro')/norm(A,'fro');
fprintf('%5.2e   : fro error of build \n',err2);

% test fast direct solver
opts.ifproxy = false;
F = chunkerflam(chnkr,K,1.0,opts);

start = tic; sol2 = rskelf_sv(F,rhs); t1 = toc(start);

if(l2scale)
    wts_rep = repmat(wts(:).', K.opdims(1),1);
    wts_rep = wts_rep(:);
    sol2 = sol2./sqrt(wts_rep);
end


fprintf('%5.2e s : time for rskelf_sv \n',t1)

err = norm(sol-sol2,'fro')/norm(sol,'fro');

fprintf('difference between fast-direct and iterative %5.2e\n',err)


