
%CHUNKERMAT_TRUEPOLYGONTEST
%
% test the matrix builder and do a basic solve for a true polygonal
% domain. This is the Neumann problem...

iseed = 8675309;
rng(iseed);


% pre-defined vertices for a barbell shape

verts = chnk.demo.barbell(2.0,2.0,1.0,1.0);

% chunk up polygon, true corner, adaptive refinement

cparams = [];
cparams.rounded = true;
cparams.depth = 20;
start = tic; chnkr = chunkerpoly(verts,cparams);
chnkr = refine(chnkr);
chnkr = sort(chnkr);
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)


wchnkr = chnkr.wts;
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

% figure(1)
% clf
% hold off
% plot(chnkr)
% hold on
% scatter(sources(1,:),sources(2,:),'o')
% scatter(targets(1,:),targets(2,:),'x')
% axis equal 

%

kerns= @(s,t) chnk.lap2d.kern(s,t,'s');
kernd= @(s,t) chnk.lap2d.kern(s,t,'d');
kernsprime = @(s,t) chnk.lap2d.kern(s,t,'sprime');

% eval du/dn on bdry

srcinfo = []; srcinfo.r = sources;
targinfo = []; targinfo.r = chnkr.r(:,:); 
targinfo.d = chnkr.d(:,:);
kernmatsprime = kernsprime(srcinfo,targinfo);
dudnbdry = kernmatsprime*strengths;

% eval u at targets

targinfo.r = targets; targinfo.d = [];
kernmatstarg = kerns(srcinfo,targinfo);
utarg = kernmatstarg*strengths;

%

% build laplace S' matrix

opts = [];
opts.l2scale = true;
start = tic; Sprime = chunkermat(chnkr,kernsprime,opts);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

wmat = onesmat(chnkr);
if opts.l2scale
    wmat = diag(sqrt(wchnkr))*wmat*diag(1./sqrt(wchnkr));
end

sys = 0.5*eye(chnkr.k*chnkr.nch) + wmat + Sprime;
%sys = 0.5*eye(chnkr.k*chnkr.nch) + Sprime;
 

rhs = dudnbdry; rhs = rhs(:);
if opts.l2scale
    rhs = rhs.*sqrt(wchnkr);
end
start = tic; sol = gmres(sys,rhs,[],1e-13,100); t1 = toc(start);

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

opts.verb=false;
start=tic; Ssol = chunkerkerneval(chnkr,kerns,sol,targets,opts); 
t1 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n',t1)

% only determined up to a constant

diffmean = mean(utarg-Ssol);
Ssol = Ssol + diffmean;

% error

relerr = norm(utarg-Ssol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
relerr2 = norm(utarg-Ssol,'inf')/dot(abs(sol(:)),wchnkr(:));
fprintf('relative frobenius error %5.2e\n',relerr);
fprintf('relative l_inf/l_1 error %5.2e\n',relerr2);

assert(relerr < 1e-10);

