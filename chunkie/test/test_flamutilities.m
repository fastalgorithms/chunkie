
%TEST_FLAMUTILITIES
%
% test the FLAM matrix builder and do a basic solve

iseed = 8675309;
rng(iseed);

addpaths_loc();

cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 5;
pref = []; 
pref.k = 30;
narms = 3;
amp = 0.25;
start = tic; chnkr = chunkfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

wts = whts(chnkr);

fprintf('%5.2e s : time to build geo\n',t1)

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
plot(chnkr)
hold on
scatter(sources(1,:),sources(2,:),'o')
scatter(targets(1,:),targets(2,:),'x')
axis equal 

%

kerns = @(s,t,sn,tn) chnk.lap2d.kern(s,t,sn,tn,'s');

% eval u on bdry

targs = chnkr.r; targs = reshape(targs,2,chnkr.k*chnkr.nch);
targstau = taus(chnkr); 
targstau = reshape(targstau,2,chnkr.k*chnkr.nch);

kernmats = kerns(sources,targs,[],targstau);
ubdry = kernmats*strengths;

% eval u at targets

kernmatstarg = kerns(sources,targets,[],[]);
utarg = kernmatstarg*strengths;

%

% build laplace dirichlet matrix

fkern = @(s,t,stau,ttau) chnk.lap2d.kern(s,t,stau,ttau,'D');
opdims(1) = 1; opdims(2) = 1;
intparams.intorder = chnkr.k;
start = tic; D = chunkskernmat(chnkr,fkern,opdims,intparams);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = -0.5*eye(chnkr.k*chnkr.nch) + D;

% build sparse tridiag part 

start = tic; spmat = chunkskernmattd(chnkr,fkern,opdims,intparams);
t1 = toc(start);
fprintf('%5.2e s : time to build tridiag\n',t1)

spmat = spmat -0.5*speye(chnkr.k*chnkr.nch);

% % test matrix entry evaluator
% 
% start = tic; 
% sys2 = kernbyindex(1:chnkr.npt,1:chnkr.npt,chnkr,wts,fkern,opdims,spmat);
% t1 = toc(start);
% 
% fprintf('%5.2e s : time for mat entry eval on whole mat\n',t1)


xflam = chnkr.r(:,:);
matfun = @(i,j) kernbyindex(i,j,chnkr,wts,fkern,opdims,spmat);
[pr,ptau,pw,pin] = proxy_square_pts();
ifaddtrans = true;
pxyfun = @(x,slf,nbr,l,ctr) proxyfun(slf,nbr,l,ctr,chnkr,wts, ...
    fkern,opdims,pr,ptau,pw,pin,ifaddtrans);
start = tic; F = rskelf(matfun,xflam,200,1e-14,pxyfun); t1 = toc(start);

afun = @(x) rskelf_mv(F,x);

fprintf('%5.2e s : time for flam compress\n',t1)

err2 = norm(sys2-sys,'fro')/norm(sys,'fro');
fprintf('%5.2e   : fro error of build \n',err2)

rhs = ubdry; rhs = rhs(:);
start = tic; sol = gmres(sys,rhs,[],1e-14,100); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

rhs = ubdry; rhs = rhs(:);
start = tic; sol3 = rskelf_sv(F,rhs); t1 = toc(start);

fprintf('%5.2e s : time for rskelf_sv \n',t1)

err = norm(sol-sol3,'fro')/norm(sol,'fro');

fprintf('difference between fast-direct and iterative %5.2e\n',err)

% evaluate at targets and compare

opts.usesmooth=false;
opts.verb=false;
opts.quadkgparams = {'RelTol',1e-16,'AbsTol',1.0e-16};
start=tic; Dsol = chunkerintkern(chnkr,fkern,opdims,sol2,targets,opts); 
t1 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n',t1)

%

wchnkr = whts(chnkr);

relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
relerr2 = norm(utarg-Dsol,'inf')/dot(abs(sol(:)),wchnkr(:));
fprintf('relative frobenius error %5.2e\n',relerr);
fprintf('relative l_inf/l_1 error %5.2e\n',relerr2);

