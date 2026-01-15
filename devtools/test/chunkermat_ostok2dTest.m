% chunkermat_ostok2dTest0();

% function chunkermat_ostok2dTest0()

%CHUNKERMAT_OSTOK2DTEST
%
% test the matrix builder and do a basic solve

% % % % % iseed = 8675309;
% % % % % rng(iseed);
% % % % % 
% % % % % cparams = [];
% % % % % cparams.eps = 1.0e-10;
% % % % % cparams.nover = 1;
% % % % % pref = []; 
% % % % % pref.k = 16;
% % % % % narms = 3;
% % % % % amp = 0.25;
start = tic;
rad = 1; ctr = [0.0;0.0];
circfun = @(t) ctr + rad*[cos(t(:).');sin(t(:).')];


opts = [];
opts.maxchunklen = 0.5;
pref = []; pref.nchmax = 100000;
chnkr1 = chunkerfunc(circfun, opts, pref);


chnkr2 = 0.2*chnkr1;
chnkr2 = chnkr2.reverse;
chnkr = merge([chnkr1, chnkr2]); % only one hole

ns = 10;
nt = 100;

sources = 2.05*circfun(1:ns);  % only one hole


targets = 0.90*circfun(1:nt);


strengths = randn(2*ns,1); % only one hole



sources_n = 0.0*ones(2,ns)+0.5*rand(2,ns);  % only one hole

plot(chnkr, 'r.'); hold on;
plot(targets(1,:), targets(2,:), 'kx');
plot(sources(1,:), sources(2,:), 'bo');hold on;

t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

% section %% 
zk = 1.3;
kernsp = kernel('ostok', 'sp', zk);
kerns = kernel('ostok', 's', zk);
kernd = kernel('ostok', 'd', zk);
kernsink = kernel('ostok', 'sinkv', zk);

% eval u on bdry


srcinfo = []; 
srcinfo.r = sources; 
srcinfo.n = sources_n;
% kernmats = kernsink.eval(srcinfo, chnkr);
% kernmats = kernsp.eval(srcinfo, chnkr); % for sprime
% kernmats = kerns.eval(srcinfo, chnkr); % for single layer
kernmats = kernd.eval(srcinfo, chnkr); % for double layer 
ubdry = kernmats*strengths; 



% section %%
ubdry2 = kernd.eval(srcinfo, chnkr)*strengths;
n1 = chnkr1.npt;
ux1 = ubdry2(1:2:2*n1);
uy1 = ubdry2(2:2:2*n1);
nx1 = chnkr1.n(1,:).';
ny1 = chnkr1.n(2,:).';

r = sum((ux1.*nx1 + uy1.*ny1).*chnkr1.wts(:))


% section %%
% eval u at targets

targinfo = []; 
targinfo.r = targets;
targets_n = rand(2, nt); 
targets_n = targets_n./sqrt(targets_n(1,:).^2+targets_n(2,:).^2);
targinfo.n = targets_n;
% kernmatstarg = kerns.eval(srcinfo, targinfo); % for sprime and single layer
kernmatstarg = kernd.eval(srcinfo, targinfo); % for double layer 
utarg = kernmatstarg*strengths; 
%%%%% exact solution uisng analyticity not code
% section %%
% solve

% build Oscillatory Stokes dirichlet matrix
% fkern = kernel('ostok', 'sp', zk); % for sprime
fkern = kernel('ostok', 'd', zk);  % for double layer
% fkern = kernel('ostok', 's', zk);  % for single layer
uu = fkern.eval(chnkr, targinfo);
% section %%
start = tic; 
D = chunkermat(chnkr, fkern);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

% sys = 0.5*eye(size(D,1)) + D;
% sys = sys + normonesmat(chnkr)/sum(chnkr.wts(:)); % for sprime

sys = -0.5*eye(size(D,1)) + D;
sys = sys + normonesmat(chnkr)/sum(chnkr.wts(:)); % for double layer


% sys = D; % for single layer

rhs = ubdry; 
rhs = rhs(:);

start = tic; 
sol = gmres(sys,rhs,[],1e-12,1000); 
t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)



% evaluate at targets and compare

opts.usesmooth=false;
opts.verb=false;

%%%%%%%% SINGLE LAYER TEST 

% fkern =  kernel('ostok', 's', zk);
% Dsol = chunkerkerneval(chnkr, fkern, sol, targets, opts);
% relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
% fprintf('relative frobenius error %5.2e\n',relerr);
% 
% fprintf('relative frobenius error %5.2e\n',relerr);

%%%%%%%% DOUBLE LAYER TEST

fkernd =  kernel('ostok', 'd', zk);
Dsol = chunkerkerneval(chnkr, fkernd, sol, targets, opts);
relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
fprintf('relative frobenius error %5.2e\n',relerr);

fprintf('relative frobenius error %5.2e\n',relerr);

%%%%%%%   SPRIME TEST

% fkerns =  kernel('ostok', 'sp', zk); % for sprime
% Dsol = chunkerkerneval(chnkr, fkerns, sol, targets, opts);
% relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
% fprintf('relative frobenius error %5.2e\n',relerr);
% 
% fprintf('relative frobenius error %5.2e\n',relerr);
% 
% assert(relerr < 1e-10);

%%%%%%%%%% SINK VELOCITY 

% end