
%FLAMUTILITIESTEST
%
% test the FLAM matrix builder and do a basic solve

clearvars; close all;
iseed = 8675309;
rng(iseed);

addpaths_loc();

nverts = 3; 
verts = exp(1i*2*pi*(0:(nverts-1))/nverts);
verts = [real(verts);imag(verts)];


iind = 1:nverts;
jind = 1:nverts;

iind = [iind iind];
jind = [jind jind + 1];
jind(jind>nverts) = 1;
svals = [-ones(1,nverts) ones(1,nverts)];
edge2verts = sparse(iind,jind,svals,nverts,nverts);

amp = 0.1;
frq = 2;
fchnks    = cell(1,size(edge2verts,1));
for icurve = 1:size(edge2verts,1)
    fchnks{icurve} = @(t) sinearc(t,amp,frq);
end
cparams = [];
cparams.nover = 2;
[cgrph] = chunkgraphinit(verts,edge2verts,fchnks,cparams);

cgrph = balance(cgrph);


wts = weights(cgrph); wts = wts(:);

% sources

ns = 10;
ts = 0.0+2*pi*rand(ns,1);
sources = 3.0*[cos(ts)';sin(ts)'];
strengths = randn(ns,1);

% targets

nt = 1;
ts = 0.0+2*pi*rand(nt,1);
targets = 0.2*[cos(ts)'; sin(ts)'];
targets = targets.*repmat(rand(1,nt),2,1);

% plot geo and sources

figure(1)
clf
hold off
plot(cgrph)
hold on
scatter(sources(1,:),sources(2,:),'o')
scatter(targets(1,:),targets(2,:),'x')
axis equal 

%

kerns = kernel('lap','s');

% eval u on bdry

srcinfo = []; srcinfo.r = sources;
targinfo = []; targinfo.r = cgrph.r(:,:); 
targinfo.d = cgrph.d(:,:);
ubdry = kerns.fmm(1e-12,srcinfo,targinfo,strengths);

% eval u at targets

targinfo = []; targinfo.r = targets;
utarg = kerns.fmm(1e-12,srcinfo,targinfo,strengths);

%

% build laplace dirichlet matrix * (-2)

kernd = -2*kernel('lap','d');

opts = [];
opts.rcip = true;
start = tic; D = chunkermat(cgrph,kernd,opts);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

npt = cgrph.npt;
sys = eye(npt) + D;

% build sparse tridiag part 

opts.nonsmoothonly = true;
opts.rcip = true;
start = tic; spmat = chunkermat(cgrph,kernd,opts);
t1 = toc(start);
fprintf('%5.2e s : time to build tridiag\n',t1)

spmat = spmat + speye(npt);

% test matrix entry evaluator
start = tic; 
% opdims = [1 1];
opdims = ones([2,1,1]);

sys2 = chnk.flam.kernbyindex(1:npt,1:npt,cgrph,kernd,opdims,spmat);



t1 = toc(start);

fprintf('%5.2e s : time for mat entry eval on whole mat\n',t1)

err2 = norm(sys2-sys,'fro')/norm(sys,'fro');
fprintf('%5.2e   : fro error of build \n',err2)


xflam = cgrph.r(:,:);
matfun = @(i,j) chnk.flam.kernbyindex(i,j,cgrph,kernd,opdims,spmat);
[pr,ptau,pw,pin] = chnk.flam.proxy_square_pts();
ifaddtrans = true;
pxyfun = @(x,slf,nbr,l,ctr) chnk.flam.proxyfun(slf,nbr,l,ctr,cgrph, ...
    kernd,opdims,pr,ptau,pw,pin,ifaddtrans);


start = tic; F = rskelf(matfun,xflam,200,1e-14,[]); t1 = toc(start);
% F = chunkerflam(cgrph,kernd,1.0);

fprintf('%5.2e s : time for flam rskelf compress\n',t1)

pxyfunr = @(rc,rx,cx,slf,nbr,l,ctr) chnk.flam.proxyfunr(rc,rx,slf,nbr,l, ...
        ctr,cgrph,kernd,opdims,pr,ptau,pw,pin);

opts = [];
start = tic; F2 = rskel(matfun,xflam,xflam,200,1e-14,pxyfunr,opts); t1 = toc(start);

fprintf('%5.2e s : time for flam rskel compress\n',t1)

afun = @(x) rskelf_mv(F,x);



rhs = ubdry; rhs = rhs(:);
start = tic; sol = gmres(sys,rhs,[],1e-14,100); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

rhs = ubdry; rhs = rhs(:);
start = tic; sol3 = rskelf_sv(F,rhs); t1 = toc(start);

fprintf('%5.2e s : time for rskelf_sv \n',t1)

err = norm(sol-sol3,'fro')/norm(sol,'fro');

fprintf('difference between fast-direct and iterative %5.2e\n',err)

% assert(err < 1e-10);

% evaluate at targets and compare

opts.usesmooth=false;
opts.verb=false;
opts.quadkgparams = {'RelTol',1e-16,'AbsTol',1.0e-16};

% Collapse cgrph into chnkrtotal
chnkrs = cgrph.echnks;
chnkrtotal = merge(chnkrs);
start=tic; Dsol = chunkerkerneval(chnkrtotal,kernd,sol,targets,opts); 
t1 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n',t1)

%

wcgrph = weights(cgrph);

relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkrtotal.nch)*norm(utarg,'fro'));
relerr2 = norm(utarg-Dsol,'inf')/dot(abs(sol(:)),wcgrph(:));
fprintf('relative frobenius error %5.2e\n',relerr);
fprintf('relative l_inf/l_1 error %5.2e\n',relerr2);



function [r,d,d2] = sinearc(t,amp,frq)
xs = t;
ys = amp*sin(frq*t);
xp = ones(size(t));
yp = amp*frq*cos(frq*t);
xpp = zeros(size(t));
ypp = -frq*frq*amp*sin(t);

r = [(xs(:)).'; (ys(:)).'];
d = [(xp(:)).'; (yp(:)).'];
d2 = [(xpp(:)).'; (ypp(:)).'];
end
