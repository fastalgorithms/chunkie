%ELASTIC_TEST
% test the kernel definitions for elasticity

clear all
iseed = 1234;
rng(iseed,'twister');
clc
addpaths_loc();

% pde parameters 

lam = 1.5;
mu = 2.1;

% set up charge outside and target inside starfish curve
% for Green's ID test

cparams = []; cparams.eps=1e-10;
chnkr = chunkerfunc(@(t) starfish(t),cparams);
chnkr = chnkr.sort();

f = [-1.3;2]; % charge strength
src = []; src.r = [2;1]; src.n = randn(2,1);
targ = []; targ.r = 0.1*randn(2,4); targ.n = randn(2,4);

%figure(1)
%clf
%plot(chnkr,'g')
%hold on
%plot(src.r(1),src.r(2),'bo')
%plot(targ.r(1),targ.r(2),'rx')
%axis equal

% test that the kernels do what we expect 
% (Green's ID, solve pde, etc)

niter=5;
[pde_errs,pdedalt_errs,div_errs,trac_errs,gid_err,daltgrad_errs] = ...
    test_kernels(chnkr,src,targ,lam,mu,f,niter);

assert(gid_err < 1e-14);
assert(min(pde_errs) < 1e-5);
assert(min(pdedalt_errs) < 1e-5);
assert(min(div_errs) < 1e-7);
assert(min(trac_errs) < 1e-7);
assert(min(daltgrad_errs) < 1e-7);

%

diag = 0.5*((lam+3*mu)/(lam+2*mu));

fkern = @(s,t) chnk.elast2d.kern(lam,mu,s,t,'dalt');
start = tic(); sysmat = chunkermat(chnkr,fkern); toc(start)
sysmat2 = sysmat + diag*eye(2*chnkr.npt);

t = []; t.r = chnkr.r(:,:); 
t.n = chnkr.n(:,:); t.d = chnkr.d(:,:);

wts = chnkr.wts; wts = wts(:); 
wts2 = [wts(:).'; wts(:).']; wts2 = wts2(:);

ubdry = elasticlet(lam,mu,src,t,f);

sigma = sysmat2\ubdry;

dmatalt = chnk.elast2d.kern(lam,mu,t,targ,'dalt');
utargsol = dmatalt*(wts2.*sigma);

utarg = elasticlet(lam,mu,src,targ,f);

err_dir = norm(utarg-utargsol)/norm(utarg);

assert(err_dir < 1e-10);

%

sp = -0.5*(mu/(lam+2*mu)); stoktrac = -0.5*(lam+mu)/(lam+2*mu);
diag = sp+stoktrac;
fkern = kernel();
fkern.eval = @(s,t) chnk.elast2d.kern(lam,mu,s,t,'strac');
fkern.sing = 'pv';
fkern.opdims = [2,2];

start = tic(); sysmat = chunkermat(chnkr,fkern); toc(start)
sysmat2 = sysmat + diag*eye(2*chnkr.npt);

t = []; t.r = chnkr.r(:,:); 
t.n = chnkr.n(:,:); t.d = chnkr.d(:,:);

wts = weights(chnkr); wts = wts(:); 
wts2 = [wts(:).'; wts(:).']; wts2 = wts2(:);

[~,tracbdry] = elasticlet(lam,mu,src,t,f);

sigma = gmres(sysmat2,tracbdry,[],1e-13,100);

smat = chnk.elast2d.kern(lam,mu,t,targ,'s');
utargsol = smat*(wts2.*sigma);

utarg = elasticlet(lam,mu,src,targ,f);

targperp = chnk.perp(targ.r);
nt = size(targ.r,2);
nspace = [repmat(eye(2),nt,1) targperp(:)];
cfs = nspace\(utarg-utargsol);
utargsol = utargsol + nspace*cfs;

err_neu = norm(utarg-utargsol)/norm(utarg);

assert(err_neu < 1e-10);

function [pde_errs,pdedalt_errs,div_errs,trac_errs,gid_err,daltgrad_errs] = ...
    test_kernels(chnkr,s,targ,lam,mu,f,niter)
%TEST_KERNELS
%
% test the elasticity kernel definitions against finite
% difference approximations of pde, divergence, traction


% check Green's ID

t = []; t.r = chnkr.r(:,:); 
t.n = chnkr.n(:,:); t.d = chnkr.d(:,:);

wts = chnkr.wts; wts = wts(:); 
wts2 = [wts(:).'; wts(:).']; wts2 = wts2(:);

[u,trac] = elasticlet(lam,mu,s,t,f);
[utarg,tractarg] = elasticlet(lam,mu,s,targ,f);
smat = chnk.elast2d.kern(lam,mu,t,targ,'s');
dmat = chnk.elast2d.kern(lam,mu,t,targ,'d');

uint = smat*(wts2.*trac) - dmat*(wts2.*u); uint = reshape(uint,size(utarg));
gid_err = norm(-ones(size(utarg))-uint./utarg,'fro');

%fprintf('%5.2e Greens ID err\n',gid_err);

% finite diff tests 

pde_errs = zeros(niter,1);
pdedalt_errs = zeros(niter,1);
div_errs = zeros(niter,1);
trac_errs = zeros(niter,1);
daltgrad_errs = zeros(niter,1);
is = randi(chnkr.npt);
s = [];
s.r = chnkr.r(:,is);
s.d = chnkr.d(:,is);
s.n = chnkr.n(:,is);
it = randi(chnkr.npt);
t = [];
t.r = chnkr.r(:,it);
t.d = chnkr.d(:,it);
t.n = chnkr.n(:,it);

[u00,trac00,dub00,dalt00,dalt00grad,dalt00trac,u00grad] = elasticlet(lam,mu,s,t,f);

for i = 1:niter
    h = 0.1^i;
    t01 = t; t10 = t; t0m1 = t; tm10=t;
    t11 = t; t1m1 = t; tm1m1 = t; tm11=t;
    t01.r(2) = t01.r(2)+h;
    t10.r(1) = t10.r(1)+h;
    t0m1.r(2) = t0m1.r(2)-h;
    tm10.r(1) = tm10.r(1)-h;
    t11.r(:) = t11.r(:)+h;
    t1m1.r(:) = t1m1.r(:)+[h;-h];
    tm1m1.r(:) = tm1m1.r(:)-h;
    tm11.r(:) = tm11.r(:)+[-h;h];
    [u01,trac01,dub01,dalt01] = elasticlet(lam,mu,s,t01,f);
    [u10,trac10,dub10,dalt10] = elasticlet(lam,mu,s,t10,f);
    [u0m1,trac0m1,dub0m1,dalt0m1] = elasticlet(lam,mu,s,t0m1,f);
    [um10,tracm10,dubm10,daltm10] = elasticlet(lam,mu,s,tm10,f);
    [u11,trac11,dub11,dalt11] = elasticlet(lam,mu,s,t11,f);
    [u1m1,trac1m1,dub1m1,dalt1m1] = elasticlet(lam,mu,s,t1m1,f);
    [um1m1,tracm1m1,dubm1m1,daltm1m1] = elasticlet(lam,mu,s,tm1m1,f);
    [um11,tracm11,dubm11,daltm11] = elasticlet(lam,mu,s,tm11,f);
    
    lapu = (u01+u10+u0m1+um10-4*u00)/h^2;
    ux = (u10-um10)/(2*h);
    uy = (u01-u0m1)/(2*h);
    uxx = (u10+um10-2*u00)/h^2;
    uyy = (u01+u0m1-2*u00)/h^2;
    uxy = (u11-um11-u1m1+um1m1)/(4*h^2);
    uyx = uxy;

    lapdub = (dub01+dub10+dub0m1+dubm10-4*dub00)/h^2;
    dubx = (dub10-dubm10)/(2*h);
    duby = (dub01-dub0m1)/(2*h);
    dubxx = (dub10+dubm10-2*dub00)/h^2;
    dubyy = (dub01+dub0m1-2*dub00)/h^2;
    dubxy = (dub11-dubm11-dub1m1+dubm1m1)/(4*h^2);
    dubyx = dubxy;
    
    lapdalt = (dalt01+dalt10+dalt0m1+daltm10-4*dalt00)/h^2;
    daltx = (dalt10-daltm10)/(2*h);
    dalty = (dalt01-dalt0m1)/(2*h);
    daltxx = (dalt10+daltm10-2*dalt00)/h^2;
    daltyy = (dalt01+dalt0m1-2*dalt00)/h^2;
    daltxy = (dalt11-daltm11-dalt1m1+daltm1m1)/(4*h^2);
    daltyx = daltxy;
    
    tracx = (trac10-tracm10)/(2*h);
    tracy = (trac01-trac0m1)/(2*h);
    
    pdeuh = mu*lapu + (lam+mu)*[uxx(1)+uxy(2);uyx(1)+uyy(2)];
    %fprintf('%5.2e pde err\n',norm(pdeuh)/norm(u00));
    pdedubh = mu*lapdub + (lam+mu)*[dubxx(1)+dubxy(2);dubyx(1)+dubyy(2)];
    %fprintf('%5.2e pde of doublet err\n',norm(pdedubh)/norm(dub00));
    pdedalth = mu*lapdalt + (lam+mu)*[daltxx(1)+daltxy(2);daltyx(1)+daltyy(2)];
    pdedalt_errs(i) = norm(pdedalth)/norm(u00);
    %fprintf('%5.2e pde of alt doublet err\n',norm(pdedalth)/norm(dalt00));
    pde_errs(i) = norm(pdeuh)/norm(u00);
    divtrach = tracx(1)+tracy(2);
    %fprintf('%5.2e div err\n',norm(divtrach)/norm(trac00))
    div_errs(i) = norm(divtrach)/norm(trac00);
    jact = [ux uy]; epsmat = 0.5*(jact + jact.');
    ugradh = jact.'; ugradh = ugradh(:);
    fprintf('%5.2e sgrad err\n',norm(ugradh-u00grad)/norm(u00grad));
    daltjact = [daltx dalty]; daltepsmat = 0.5*(daltjact + daltjact.');
    tracuh = (lam*(ux(1)+uy(2))*eye(2) + 2*mu*epsmat)*t.n;
    tracdalt = (lam*(daltx(1)+dalty(2))*eye(2) + 2*mu*daltepsmat)*t.n;
    %fprintf('%5.2e trac err\n',norm(tracuh-trac00)/norm(trac00));
    fprintf('%5.2e dalttrac err\n',norm(tracdalt-dalt00trac)/norm(dalt00trac));
    trac_errs(i) = norm(tracuh-trac00)/norm(trac00);
    
    daltgrad_errs(i) = ...
        norm(dalt00grad(1:2:end)-daltx)/norm(dalt00grad(1:2:end)) ...
        + norm(dalt00grad(2:2:end)-dalty)/norm(dalt00grad(2:2:end));
    
end

end

function [u,trac,dub,dalt,daltgrad,dalttrac,ugrad] = elasticlet(lam,mu,s,t,f)

mat = chnk.elast2d.kern(lam,mu,s,t,'s');
mattrac = chnk.elast2d.kern(lam,mu,s,t,'strac');
matdub = chnk.elast2d.kern(lam,mu,s,t,'d');
matdalt = chnk.elast2d.kern(lam,mu,s,t,'dalt');
matdaltgrad = chnk.elast2d.kern(lam,mu,s,t,'daltgrad');
matdalttrac = chnk.elast2d.kern(lam,mu,s,t,'dalttrac');
matgrad = chnk.elast2d.kern(lam,mu,s,t,'sgrad');
u = mat*f;
trac = mattrac*f;
dub = matdub*f;
dalt = matdalt*f;
daltgrad = matdaltgrad*f;
dalttrac = matdalttrac*f;
ugrad = matgrad*f;

end

