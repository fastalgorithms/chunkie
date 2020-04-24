%DEMO_SCATTER
%
% Define a scattering problem on a starfish-shaped domain and solve
%

iseed = 8675309;
rng(iseed);
addpaths_loc();

% planewave vec

kvec = 10*[1;-1.5];

%

zk = norm(kvec);

%   

cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 0;
cparams.maxchunklen = 4.0/zk;
pref = []; 
pref.k = 16;
narms =5;
amp = 0.25;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

[~,~,info] = sortinfo(chnkr);
assert(info.ier == 0);

% plot geometry and data

figure(1)
clf
plot(chnkr,'-x')
hold on
quiver(chnkr)
axis equal

%


% solve and visualize the solution

% build CFIE

fkern = @(s,t,stau,ttau) chnk.helm2d.kern(zk,s,t,stau,ttau,'c',1);
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; sysmat = chunkermat(chnkr,fkern,opts);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = 0.5*eye(chnkr.k*chnkr.nch) + sysmat;

rhs = -planewave(kvec(:),chnkr.r(:,:)); rhs = rhs(:);
start = tic; sol = gmres(sys,rhs,[],1e-14,100); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

% evaluate at targets and plot

rmin = min(chnkr); rmax = max(chnkr);
xl = rmax(1)-rmin(1);
yl = rmax(2)-rmin(2);
nplot = 300;
xtarg = linspace(rmin(1)-xl,rmax(1)+xl); 
ytarg = linspace(rmin(2)-yl,rmax(2)+yl,nplot);
[xxtarg,yytarg] = meshgrid(xtarg,ytarg);
targets = zeros(2,length(xxtarg(:)));
targets(1,:) = xxtarg(:); targets(2,:) = yytarg(:);

start = tic;
chnkr2 = chnkr;
chnkr2 = chnkr2.makedatarows(1);
chnkr2.data(1,:) = sol(:);
optref = []; optref.nover = 4;
chnkr2 = chnkr2.refine(optref);
sol2 = chnkr2.data(1,:);
t1 = toc(start);

fprintf('%5.2e s : time to oversample boundary\n',t1)

%

start = tic; in = chunkerinteriorflam(chnkr,targets); t1 = toc(start);
out = ~in;

fprintf('%5.2e s : time to find points in domain\n',t1)

% compute layer potential based on oversample boundary

wts2 = weights(chnkr2);

matfun = @(i,j) kernbyindexr(i,j,targets(:,out),chnkr2,wts2,fkern,opdims);
[pr,ptau,pw,pin] = proxy_square_pts();

pxyfun = @(rc,rx,cx,slf,nbr,l,ctr) proxyfunr(rc,rx,slf,nbr,l,ctr,chnkr2,wts2, ...
    fkern,opdims,pr,ptau,pw,pin);

xflam = chnkr2.r(:,:);

fmmopts = []; fmmopts.store = 'A'; % store everything (faster apply)
start = tic; F = ifmm(matfun,targets(:,out),xflam,200,1e-14,pxyfun,fmmopts); 
t1 = toc(start);
fprintf('%5.2e s : time for ifmm form (for plotting)\n',t1)
start = tic;
uscat = ifmm_mv(F,sol2(:),matfun); t1 = toc(start);
fprintf('%5.2e s : time for ifmm apply (for plotting)\n',t1)

uin = planewave(kvec,targets(:,out));
utot = uscat(:)+uin(:);

%

maxin = max(abs(uin(:)));
maxsc = max(abs(uin(:)));
maxtot = max(abs(uin(:)));

maxu = max(max(maxin,maxsc),maxtot);

figure(2)
clf
subplot(1,3,1)
zztarg = nan(size(xxtarg));
zztarg(out) = uin;
h=pcolor(xxtarg,yytarg,imag(zztarg));
set(h,'EdgeColor','none')
hold on
plot(chnkr,'LineWidth',2)
axis equal
axis tight
colormap(redblue)
caxis([-maxu,maxu])
title('$u_{in}$','Interpreter','latex','FontSize',24)


subplot(1,3,2)
zztarg = nan(size(xxtarg));
zztarg(out) = uscat;
h=pcolor(xxtarg,yytarg,imag(zztarg));
set(h,'EdgeColor','none')
hold on
plot(chnkr,'LineWidth',2)
axis equal
axis tight
colormap(redblue)
caxis([-maxu,maxu])
title('$u_{scat}$','Interpreter','latex','FontSize',24)

subplot(1,3,3)
zztarg = nan(size(xxtarg));
zztarg(out) = utot;
h=pcolor(xxtarg,yytarg,imag(zztarg));
set(h,'EdgeColor','none')
hold on
plot(chnkr,'LineWidth',2)
axis equal
axis tight
colormap(redblue)
caxis([-maxu,maxu])
title('$u_{tot}$','Interpreter','latex','FontSize',24)
