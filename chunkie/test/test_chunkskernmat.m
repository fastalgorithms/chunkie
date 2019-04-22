
%TEST_CHUNKSKERNMAT
%
% test the matrix builder and do a basic solve

iseed = 8675309;
rng(iseed);

addpaths_loc();

cparams = [];
cparams.eps = 1.0e-4;
cparams.nover = 0;
pref = []; 
pref.k = 16;
narms = 10;
amp = 0.5;
start = tic; chnkr = chunkfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

figure(1)
plot(chnkr)
hold on
quiver(chnkr)
axis equal

fprintf('%5.2e s : time to build geo\n',t1)



% sources

ns = 10;
ts = 0.0+2*pi*rand(ns,1);
sources = starfish(ts,narms,amp);
sources = 3.0*sources;
strengths = randn(ns,1);

% targets

nt = 10;
ts = 0.0+2*pi*rand(nt,1);
targets = starfish(ts,narms,amp);
targets = targets.*repmat(rand(1,nt),2,1);

% plot geo and sources

xs = chnkr.r(1,:,:); xmin = min(xs(:)); xmax = max(xs(:));
ys = chnkr.r(2,:,:); ymin = min(ys(:)); ymax = max(ys(:));

hold off
plot(chnkr)
hold on
scatter(sources(1,:),sources(2,:),'o')
scatter(targets(1,:),targets(2,:),'x')
axis equal 

%

kerns = @(s,t,sn,tn) glapkern(s,t,sn,tn,'s');

% eval u on bdry

targs = chnkr.r; targs = reshape(targs,2,chnkr.k*chnkr.nch);
targstau = taus(chnkr); 
targstau = reshape(targstau,2,chnkr.k*chnkr.nch);

kernmats = kerns(sources,targs,[],targstau);
ubdry = kernmats*strengths;

% eval u at targets

kernmatstarg = kerns(sources,targets,[],[]);
utarg = kernmatstarg*strengths;

%%

% build laplace dirichlet matrix

fkern = @(s,t,stau,ttau) glapkern(s,t,stau,ttau,'D');
opdims(1) = 1; opdims(2) = 1;
intparams.intorder = chnkr.k;
start = tic; D = chunkskernmat(chnkr,fkern,opdims,intparams);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = -0.5*eye(chnkr.k*chnkr.nch) + D;

rhs = ubdry; rhs = rhs(:);
start = tic; sol = gmres(sys,rhs,[],1e-14,100); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

start = tic; sol2 = sys\rhs; t1 = toc(start);

fprintf('%5.2e s : time for dense backslash solve\n',t1)

err = norm(sol-sol2,'fro')/norm(sol2,'fro');

fprintf('difference between direct and iterative %5.2e\n',err)

% evaluate at targets and compare

opts.usesmooth=false;
opts.verb=false;
opts.quadkgparams = {'RelTol',1.0e-15,'AbsTol',1.0e-15};
start=tic; Dsol = chunkerintkern(chnkr,fkern,opdims,sol,targets,opts); 
t1 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n',t1)

%

relerr = norm(utarg-Dsol,'fro')/norm(utarg,'fro');

fprintf('relative frobenius error %5.2e\n',relerr);

