
%TEST_TRAPKERNMAT
%
% test the matrix builder and do a basic solve

zk = 1.0;
quadorder = 16;
iseed = 8675309;
rng(iseed);

addpaths_loc();

cparams = [];
cparams.npt = 200;
pref = []; 
narms = 5;
amp = 0.25;
start = tic; trap = trapfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

%

% sources

ns = 10;
ts = 0.0+2*pi*rand(ns,1);
sources = starfish(ts,narms,amp);
sources = 3.0*sources;
strengths = randn(ns,1);

% targets

nt = 1;
ts = 0.0+2*pi*rand(nt,1);
targets = starfish(ts,narms,amp);
targets = targets.*repmat(rand(1,nt),2,1);

% plot geo and sources


figure(1)
clf
hold off
plot(trap)
hold on
quiver(trap)
scatter(sources(1,:),sources(2,:),'o')
scatter(targets(1,:),targets(2,:),'x')
axis equal 

%

kerns = @(s,t,sn,tn) chnk.lap2d.kern(s,t,sn,tn,'s');
%kerns = @(s,t,sn,tn) chnk.helm2d.kern(zk,s,t,sn,tn,'s');

% eval u on bdry

targs = trap.r;
targstau = taus(trap); 

kernmats = kerns(sources,targs,[],targstau);
ubdry = kernmats*strengths;

% eval u at targets

kernmatstarg = kerns(sources,targets,[],[]);
utarg = kernmatstarg*strengths;

%

% build laplace dirichlet matrix

fkern = @(s,t,stau,ttau) chnk.lap2d.kern(s,t,stau,ttau,'D');
%fkern = @(s,t,stau,ttau) chnk.helm2d.kern(zk,s,t,stau,ttau,'D');
opdims(1) = 1; opdims(2) = 1;

%start = tic; D = chnk.quadba.buildmat(trap,fkern,quadorder,opdims,'log');
opts = [];
opts.quad = 'balog';
opts.quadorder = quadorder;
start = tic; D = trapmat(trap,fkern,opts);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = -0.5*eye(trap.npt) + D;

rhs = ubdry; rhs = rhs(:);
start = tic; sol = gmres(sys,rhs,[],1e-14,100); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

start = tic; sol2 = sys\rhs; t1 = toc(start);

fprintf('%5.2e s : time for dense backslash solve\n',t1)

err = norm(sol-sol2,'fro')/norm(sol2,'fro');

fprintf('difference between direct and iterative %5.2e\n',err)

% evaluate at targets and compare

Dsol = trapperintkern(trap,fkern,[],sol,targets,[]);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n',t1)

%

wchnkr = whts(trap);

relerr = norm(utarg-Dsol,'fro')/(sqrt(trap.npt)*norm(utarg,'fro'));
relerr2 = norm(utarg-Dsol,'inf')/dot(abs(sol(:)),wchnkr(:));
fprintf('relative frobenius error %5.2e\n',relerr);
fprintf('relative l_inf/l_1 error %5.2e\n',relerr2);

