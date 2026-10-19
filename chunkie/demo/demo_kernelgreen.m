%DEMO_KERNELGREEN
%
% Demonstrates a kernel built using kernel.green.
% Simulates scattering from a starfish for several incident waves
%

zk = 15;
ndens = 3;            % number of incident directions
th = [-pi/3 0 pi/4];
kvecs = zk*[cos(th); sin(th)];

% discretize domain

cparams = [];
cparams.eps = 1.0e-6;
cparams.nover = 0;
cparams.maxchunklen = 4.0/zk;

pref = [];
pref.k = 16;
narms = 5;
amp = 0.25;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref);
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

% Build the kernel from the Green's function and its fmm

green = @(s,t) chnk.helm2d.green(zk,s,t);
fmm = @(eps,s,pg,t,pgt) hfmm2d(eps,zk,s,pg,t,pgt);
coefs = [1,-zk*1i];
fkern = kernel.green(green,'c',coefs,'fmm',fmm);

% build system matrix and solve for all the incident waves at once

start = tic; sysmat = chunkermat(chnkr,fkern);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = 0.5*eye(chnkr.k*chnkr.nch) + sysmat;

rhs = zeros(chnkr.npt,ndens);
for j = 1:ndens
    rhs(:,j) = -reshape(planewave(kvecs(:,j),chnkr.r(:,:)),[],1);
end
start = tic; sol = sys\rhs; t1 = toc(start);

fprintf('%5.2e s : time for dense solve\n',t1)

% evaluate at targets, all the densities in one call, with the fmm

rmin = min(chnkr); rmax = max(chnkr);
xl = rmax(1)-rmin(1);
yl = rmax(2)-rmin(2);
nplot = 300;
xtarg = linspace(rmin(1)-xl,rmax(1)+xl,nplot);
ytarg = linspace(rmin(2)-yl,rmax(2)+yl,nplot);
[xxtarg,yytarg] = meshgrid(xtarg,ytarg);
targets = zeros(2,length(xxtarg(:)));
targets(1,:) = xxtarg(:); targets(2,:) = yytarg(:);

start = tic; in = chunkerinterior(chnkr,{xtarg,ytarg}); t1 = toc(start);
out = ~in;

fprintf('%5.2e s : time to find points in domain\n',t1)

start = tic;
uscat = chunkerkerneval(chnkr,fkern,sol,targets(:,out)); t1 = toc(start);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t1)

% plot the fields

figure(1)
clf
t = tiledlayout(1,ndens,'TileSpacing','compact');

for j = 1:ndens
    utot = uscat(:,j) + reshape(planewave(kvecs(:,j),targets(:,out)),[],1);
    nexttile
    zztarg = nan(size(xxtarg));
    zztarg(out) = utot;
    h = pcolor(xxtarg,yytarg,imag(zztarg));
    set(h,'EdgeColor','none')
    clim([-2,2])
    colormap(redblue);
    hold on
    plot(chnkr,'k','LineWidth',2)
    axis equal tight
    set(gca,"box","off","Xtick",[],"Ytick",[]);
    title(sprintf('$\\theta = %.2f$',th(j)), ...
        'Interpreter','latex','FontSize',12)
end
