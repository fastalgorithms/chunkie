%CHUNKERMAT_QUADADAP_CLOSETOTOUCHINGTEST
%
% - set up two discs which are nearly touching
% - compute the solution of an exterior Dirichlet problem using 
%   Laplace combined layer potential
% - compare accuracy using the adaptive routine to do nearby panel
%   evaluation vs smooth rule for everything but self and neighbors

clearvars; close all;
addpaths_loc();
iseed = 8675309;
rng(iseed);

% define geometry


circle = @(t,ctr,rad) chnk.curves.bymode(t,rad*ones(1,1),ctr);

ctr1 = [-1;0]; ctr2 = [1;0];
rad = 0.99;

cparams = [];
cparams.nover = 4;
cparams.ifclosed = 1;
chnkr1 = chunkerfunc(@(t) circle(t,ctr1,rad),cparams);
chnkr2 = chunkerfunc(@(t) circle(t,ctr2,rad),cparams);

chnkrs = [chnkr1 chnkr2];

chnkr = merge(chnkrs);

% set up sources for exterior problem

srcs = [ (0.2*randn(2,10) + ctr1), (0.2*randn(2,10) + ctr2)];
strengths = randn(size(srcs,2),1);

% targets for exterior problems

targs = [ 0.001*randn(2,1), ([-2.5;0]+0.2*randn(2,2)),([2.5;0]+0.2*randn(2,2))];

% plot geometry

% figure(1)
% clf
% plot(chnkrs,'-r')
% hold on
% %quiver(chnkrs,'r')
% scatter(srcs(1,:),srcs(2,:),'bo')
% scatter(targs(1,:),targs(2,:),'gx')
% axis equal

%

kerns = @(s,t) chnk.lap2d.kern(s,t,'s');

% eval u on bdry

srcinfo = []; srcinfo.r = srcs; targinfo = []; targinfo.r = chnkr.r(:,:);
targinfo.d = chnkr.d(:,:);
kernmats = kerns(srcinfo,targinfo);
ubdry = kernmats*strengths;

% eval u at targets

targinfo = []; targinfo.r = targs;
kernmatstarg = kerns(srcinfo,targinfo);
utarg = kernmatstarg*strengths;

%

% build laplace dirichlet matrix (combined)
% use adaptive routine to build matrix (self done by ggq, nbor by adaptive)

eta = 1;
fkern = @(s,t) chnk.lap2d.kern(s,t,'C',eta);

type = 'log';
opts = []; opts.robust = true;
opdims = [1 1];
start = tic;
mata = chnk.quadadap.buildmat(chnkr,fkern,opdims,type,opts);
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix (adaptive)\n',t1)
start = tic; mato = chunkermat(chnkr,fkern);
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix (original GGQ)\n',t1)
start = tic; 
opts.forcepquad = true; targs_bd = [chnkr.r(1,:);chnkr.r(2,:)];
matho = chunkerkernevalmat(chnkr,fkern,targs_bd,opts); 
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix (Helsing-Ojala)\n',t1)

sysa = 0.5*eye(chnkr.k*chnkr.nch) + mata;
syso = 0.5*eye(chnkr.k*chnkr.nch) + mato;
sysho = matho;

%

rhs = ubdry(:); 
start = tic; sola = sysa\rhs; t1 = toc(start);
start = tic; solo = syso\rhs; t1 = toc(start);
start = tic; solho = sysho\rhs; t1 = toc(start);

fprintf('%5.2e s : time for dense backslash solve\n',t1)

err = norm(sola-solo,'fro')/norm(sola,'fro');

fprintf('difference between adaptive and non-adaptive %5.2e\n',err)

err = norm(sola-solho,'fro')/norm(sola,'fro');

fprintf('difference between adaptive and Helsing-Ojala %5.2e\n',err)

% evaluate at targets and compare

opts.usesmooth=false;
opts.verb=false;
opts.forcepquad = false;
opts.quadkgparams = {'RelTol',1e-16,'AbsTol',1.0e-16};
start=tic; layersola = chunkerkerneval(chnkr,fkern,sola,targs,opts); 
t1 = toc(start);
start=tic; layersolo = chunkerkerneval(chnkr,fkern,solo,targs,opts); 
t1 = toc(start);
opts.forcepquad = true;
start=tic; layersolho = chunkerkerneval(chnkr,fkern,solho,targs,opts); 
t1 = toc(start);
%

wchnkr = weights(chnkr);

relerr = norm(utarg-layersola,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
relerr2 = norm(utarg-layersola,'inf')/dot(abs(sola(:)),wchnkr(:));
fprintf('relative frobenius error %5.2e (adap)\n',relerr);
fprintf('relative l_inf/l_1 error %5.2e (adap)\n',relerr2);

assert(relerr < 1e-10);

relerr = norm(utarg-layersolo,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
relerr2 = norm(utarg-layersolo,'inf')/dot(abs(solo(:)),wchnkr(:));
fprintf('relative frobenius error %5.2e (no adap)\n',relerr);
fprintf('relative l_inf/l_1 error %5.2e (no adap)\n',relerr2);

relerr = norm(utarg-layersolho,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
relerr2 = norm(utarg-layersolho,'inf')/dot(abs(solho(:)),wchnkr(:));
fprintf('relative frobenius error %5.2e (Helsing-Ojala)\n',relerr);
fprintf('relative l_inf/l_1 error %5.2e (Helsing-Ojala)\n',relerr2);

% plot density
% 
% chnkr = chnkr.makedatarows(1);
% chnkr.data(1,:) = sola;
% 
% figure(2)
% plot3(chnkr,1)
% 
