
%TEST_CHUNKERMAT_TRUE_POLYGON
%
% test the matrix builder and do a basic solve for a true polygonal
% domain. This is the Neumann problem...

iseed = 8675309;
rng(iseed);

addpaths_loc();

% pre-defined vertices for a barbell shape

verts = chnk.demo.barbell(2.0,2.0,1.0,1.0);

% chunk up polygon, true corner, adaptive refinement

cparams = [];
cparams.rounded = false;
cparams.depth = 10;
start = tic; chnkr = chunkpoly(verts,cparams);
chnkr = refine(chnkr);
chnkr = sort(chnkr);
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

wchnkr = whts(chnkr);
wchnkr = wchnkr(:);

% sources

ns = 10;
sources = -1 + 2*(rand(2,ns) > 0.5);
sources = 3.0*sources;
strengths = randn(ns,1);

% targets

nt = 3;
targets = zeros(2,nt);
targets(:,1) = [0;0] + randn(2,1)*0.2;
targets(:,2) = [1.5;0.25] + randn(2,1)*0.2;
targets(:,3) = [-1.5;-0.25] + randn(2,1)*0.2;

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

%

kerns= @(s,t,sn,tn) chnk.lap2d.kern(s,t,sn,tn,'s');
kernsprime = @(s,t,sn,tn) chnk.lap2d.kern(s,t,sn,tn,'sprime');

% eval du/dn on bdry

targs = chnkr.r; targs = reshape(targs,2,chnkr.k*chnkr.nch);
targstau = taus(chnkr); 
targstau = reshape(targstau,2,chnkr.k*chnkr.nch);

kernmatsprime = kernsprime(sources,targs,[],targstau);
dudnbdry = kernmatsprime*strengths;

% eval u at targets

kernmatstarg = kerns(sources,targets,[],[]);
utarg = kernmatstarg*strengths;

%

% build laplace S' matrix

opts = [];
opts.l2scale = true;
start = tic; Sprime = chunkermat(chnkr,kernsprime,opts);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = 0.5*eye(chnkr.k*chnkr.nch) + Sprime;

rhs = dudnbdry; rhs = rhs(:);
if opts.l2scale
    rhs = rhs.*sqrt(wchnkr);
end
start = tic; sol = gmres(sys,rhs,[],1e-14,100); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

start = tic; sol2 = sys\rhs; t1 = toc(start);

fprintf('%5.2e s : time for dense backslash solve\n',t1)

err = norm(sol-sol2,'fro')/norm(sol2,'fro');

fprintf('difference between direct and iterative %5.2e\n',err)

if opts.l2scale % unscale solution
    sol2 = sol2./sqrt(wchnkr);
    sol = sol./sqrt(wchnkr);      
end

% evaluate at targets and compare

opts.usesmooth=false;
opts.verb=false;
opts.quadkgparams = {'RelTol',1e-16,'AbsTol',1.0e-16};
start=tic; Ssol = chunkerkerneval(chnkr,kerns,sol2,targets,opts); 
t1 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n',t1)

%


relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
relerr2 = norm(utarg-Dsol,'inf')/dot(abs(sol(:)),wchnkr(:));
fprintf('relative frobenius error %5.2e\n',relerr);
fprintf('relative l_inf/l_1 error %5.2e\n',relerr2);

