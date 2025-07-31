%DEMO_FREE_PLATE_SCATTER
%
% Define an interior scattering problem on a starfish-shaped domain and 
% solve
%

zk = 10;
nu = 0.3;

% discretize domain

cparams = [];
cparams.eps = 1.0e-8;
cparams.nover = 0;
cparams.maxchunklen = 4.0/zk; % setting a chunk length helps when the
                              % frequency is known
                              
pref = []; 
pref.k = 16;
narms1 = 5;
narms2 = 3;
amp = 0.25;
start = tic; 
chnkr1 = chunkerfunc(@(t) starfish(t,narms1,amp),cparams,pref); 
chnkr2 = chunkerfunc(@(t) starfish(t,narms2,amp,[],[],0.5),cparams,pref); 
chnkr1 = reverse(chnkr1);
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

% plot geometry and data

figure(1)
clf
plot(chnkr1,'-')
plot(chnkr2,'-')
hold on
quiver(chnkr1)
quiver(chnkr2)
axis equal

%%
chnkr = chnkr1;
start = tic;
% [sys,H] = free_sysmat(chnkr,zk,nu);
sys = clamped_sysmat(chnkr,zk);
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)
%%
% building RHS
src =[]; src.r = [2;0];
% src =[]; src.r = [0.5;0];
free_bc_kern = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_bcs', nu);
clamp_bc_kern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_bcs', nu);

flex_free_kern = @(s,t) chnk.flex2d.kern(zk, s, t, 's');

% rhs = free_bc_kern(src,chnkr);
rhs = clamp_bc_kern(src,chnkr);

% Solving linear system

% start = tic; sol = gmres(sys,rhs,[],1e-12,500); t1 = toc(start);
% fprintf('%5.2e s : time for dense gmres\n',t1)    
sol = sys\rhs;

% evaluate at targets and plot

rmin = min(chnkr); rmax = max(chnkr);
xl = rmax(1)-rmin(1);
yl = rmax(2)-rmin(2);
nplot = 100;
xtarg = linspace(rmin(1)-xl/10,rmax(1)+xl/10,nplot); 
ytarg = linspace(rmin(2)-yl/10,rmax(2)+yl/10,nplot);
[xxtarg,yytarg] = meshgrid(xtarg,ytarg);
targets = zeros(2,length(xxtarg(:)));
targets(1,:) = xxtarg(:); targets(2,:) = yytarg(:);

start = tic; in = ~chunkerinterior(reverse(chnkr),{xtarg,ytarg}); t1 = toc(start);
out = ~in;

fprintf('%5.2e s : time to find points in domain\n',t1)

start1 = tic;
% uscat = free_eval(chnkr,targets(:,out),sol,H,zk,nu);
uscat = clamp_eval(chnkr,targets(:,out),sol,zk);
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

%%
uin = flex_free_kern(src,struct("r",targets(:,out)));
utot = uscat(:)-uin(:);

maxu = max(abs(uin(:)));
maxu = max(abs(uscat(:)));

figure(2)
clf

t = tiledlayout(1,2,'TileSpacing','compact');

nexttile
zztarg = nan(size(xxtarg))+NaN*1i;
zztarg(out) = uin;
h=pcolor(xxtarg,yytarg,imag(zztarg)); h.FaceColor="interp";
set(h,'EdgeColor','none')
clim([-maxu,maxu])
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('$u^{\textrm{inc}}$','Interpreter','latex','FontSize',12)

nexttile
zztarg = nan(size(xxtarg))+NaN*1i;
zztarg(out) = uscat;
h=pcolor(xxtarg,yytarg,imag(zztarg)); h.FaceColor="interp";
set(h,'EdgeColor','none')
clim([-maxu,maxu])
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('$u^{\textrm{scat}}$','Interpreter','latex','FontSize',12)

colorbar

figure(4)
clf
t = tiledlayout(1,3,'TileSpacing','compact');

nexttile
zztarg = nan(size(xxtarg))+NaN*1i;
zztarg(out) = uin;
h=pcolor(xxtarg,yytarg,real(zztarg)); h.FaceColor="interp";
set(h,'EdgeColor','none')
clim([-maxu,maxu])
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('$u^{\textrm{inc}}$','Interpreter','latex','FontSize',12)

nexttile
zztarg = nan(size(xxtarg))+NaN*1i;
zztarg(out) = uscat;
h=pcolor(xxtarg,yytarg,real(zztarg)); h.FaceColor="interp";
set(h,'EdgeColor','none')
clim([-maxu,maxu])
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('$u^{\textrm{scat}}$','Interpreter','latex','FontSize',12)

nexttile
zztarg = nan(size(xxtarg))+NaN*1i;
zztarg(out) = utot;
h=pcolor(xxtarg,yytarg,real(zztarg)); h.FaceColor="interp";
set(h,'EdgeColor','none')
clim([-maxu,maxu])
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('$u^{\textrm{tot}}$','Interpreter','latex','FontSize',12)

colorbar

figure(3);clf
zztarg = nan(size(xxtarg))+NaN*1i;
zztarg(out) = utot;
h=pcolor(xxtarg,yytarg,log10(abs(zztarg))); h.FaceColor="interp";
set(h,'EdgeColor','none')
% clim([-maxu,maxu])
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('$\log_{10}$ error','Interpreter','latex','FontSize',12)
title(t,"Free Plate BCs")
colorbar


function [sys,H] = free_sysmat(chnkr,zk,nu)

% defining free plate kernels

fkern1 =  @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate', nu);        % build the desired kernel
double = @(s,t) chnk.lap2d.kern(s,t,'d');
hilbert = @(s,t) chnk.lap2d.kern(s,t,'hilb');

opts = [];
opts.sing = 'log';

opts2 = [];
opts2.sing = 'pv';

% building system matrix

sysmat1 = chunkermat(chnkr,fkern1, opts);
D = chunkermat(chnkr, double, opts);
H = chunkermat(chnkr, hilbert, opts2);     

sysmat = zeros(2*chnkr.npt);
sysmat(1:2:end,1:2:end) = sysmat1(1:4:end,1:2:end) + sysmat1(3:4:end,1:2:end)*H  - 2*((1+nu)/2)^2*D*D;
sysmat(2:2:end,1:2:end) = sysmat1(2:4:end,1:2:end) + sysmat1(4:4:end,1:2:end)*H;
sysmat(1:2:end,2:2:end) = sysmat1(1:4:end,2:2:end) + sysmat1(3:4:end,2:2:end);
sysmat(2:2:end,2:2:end) = sysmat1(2:4:end,2:2:end) + sysmat1(4:4:end,2:2:end);

D = [-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2];  % jump matrix 
D = kron(eye(chnkr.npt), D);

sys =  D + sysmat;

end

function sys = clamped_sysmat(chnkr,zk)
% assembling system matrix

fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate');

kappa = signed_curvature(chnkr);
kappa = kappa(:);

opts = [];
opts.sing = 'log';

start = tic;
sys = chunkermat(chnkr,fkern, opts);
sys = sys - 0.5*eye(2*chnkr.npt);
sys(2:2:end,1:2:end) = sys(2:2:end,1:2:end) + kappa.*eye(chnkr.npt);

t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)
end

function uscat = free_eval(chnkr,targets,sol,H,zk,nu)

ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_eval', nu); 

dens_comb = zeros(3*chnkr.npt,1);
dens_comb(1:3:end) = sol(1:2:end);
dens_comb(2:3:end) = H*sol(1:2:end);
dens_comb(3:3:end) = sol(2:2:end);

uscat = chunkerkerneval(chnkr, ikern,dens_comb,targets);
end

function uscat = clamp_eval(chnkr,targets,sol,zk)

ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');   

uscat = chunkerkerneval(chnkr, ikern,sol,targets);
end