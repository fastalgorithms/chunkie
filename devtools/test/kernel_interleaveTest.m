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
start = tic; chnkr = chunkerfunc(@(t) starfish(t, narms, amp), cparams, pref); 
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

zk = 1.1;

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
rhs = zeros(3*npts, 1);
rhs(1:3:end) = ubdry;

% Form matrix
A = chunkermat(chnkr, K) + eye(3*npts);
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
start = tic;
Dsol = chunkerkerneval(chnkr, Keval, sol, targets, opts);
t2 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n', t2)

wchnkr = weights(chnkr);
wchnkr = repmat(wchnkr(:).', 3, 1);
relerr  = norm(utarg-Dsol) / (sqrt(chnkr.nch)*norm(utarg));
relerr2 = norm(utarg-Dsol, 'inf') / dot(abs(sol(:)), wchnkr(:));
fprintf('relative frobenius error %5.2e\n', relerr);
fprintf('relative l_inf/l_1 error %5.2e\n', relerr2);
