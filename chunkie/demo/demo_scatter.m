%DEMO_SCATTER
%
% Define an exterior scattering problem on a starfish-shaped domain and 
% solve
%

clearvars; close all;
iseed = 8675309;
rng(iseed);
addpaths_loc();

% planewave vec

kvec = 10*[1;-1.5];

%

zk = norm(kvec);

% discretize domain

cparams = [];
cparams.eps = 1.0e-6;
cparams.nover = 0;
cparams.maxchunklen = 4.0/zk; % setting a chunk length helps when the
                              % frequency is known
                              
pref = []; 
pref.k = 16;
narms =5;
amp = 0.25;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

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

fkern = kernel('helm','c',zk,[1,-zk*1i]);
start = tic; sysmat = chunkermat(chnkr,fkern);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = 0.5*eye(chnkr.k*chnkr.nch) + sysmat;

rhs = -planewave(kvec(:),chnkr.r(:,:)); rhs = rhs(:);
start = tic; sol = gmres(sys,rhs,[],1e-13,100); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

% evaluate at targets and plot

rmin = min(chnkr); rmax = max(chnkr);
xl = rmax(1)-rmin(1);
yl = rmax(2)-rmin(2);
nplot = 400;
xtarg = linspace(rmin(1)-xl,rmax(1)+xl,nplot); 
ytarg = linspace(rmin(2)-yl,rmax(2)+yl,nplot);
[xxtarg,yytarg] = meshgrid(xtarg,ytarg);
targets = zeros(2,length(xxtarg(:)));
targets(1,:) = xxtarg(:); targets(2,:) = yytarg(:);

%

start = tic; in = chunkerinterior(chnkr,{xtarg,ytarg}); t1 = toc(start);
out = ~in;

fprintf('%5.2e s : time to find points in domain\n',t1)

% compute layer potential based on oversample boundary

start = tic;
uscat = chunkerkerneval(chnkr,fkern,sol,targets(:,out)); t1 = toc(start);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t1)

uin = planewave(kvec,targets(:,out));
utot = uscat(:)+uin(:);

%

maxin = max(abs(uin(:)));
maxsc = max(abs(uin(:)));
maxtot = max(abs(uin(:)));

maxu = max(max(maxin,maxsc),maxtot);

figure(2)
clf

t = tiledlayout(1,3,'TileSpacing','compact');

nexttile
zztarg = nan(size(xxtarg));
zztarg(out) = uin;
h=pcolor(xxtarg,yytarg,imag(zztarg));
set(h,'EdgeColor','none')
clim([-maxu,maxu])
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('$u^{\textrm{inc}}$','Interpreter','latex','FontSize',12)

nexttile
zztarg = nan(size(xxtarg));
zztarg(out) = uscat;
h=pcolor(xxtarg,yytarg,imag(zztarg));
set(h,'EdgeColor','none')
clim([-maxu,maxu])
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);

title('$u^{\textrm{scat}}$','Interpreter','latex','FontSize',12)

nexttile
zztarg = nan(size(xxtarg));
zztarg(out) = utot;
h=pcolor(xxtarg,yytarg,imag(zztarg));
set(h,'EdgeColor','none')
clim([-maxu,maxu])
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);

title('$u^{\textrm{tot}}$','Interpreter','latex','FontSize',12)
