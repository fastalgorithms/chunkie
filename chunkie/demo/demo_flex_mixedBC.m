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
narms1 = 4;
narms2 = 3;
amp = 0.25*0.5;
start = tic; 
chnkr_f = chunkerfunc(@(t) starfish(t,narms1,amp,[0;4],[],1),cparams,pref); 
chnkr_c = chunkerfunc(@(t) starfish(t,narms2,amp),cparams,pref); 

% chnkr_f = chunkerfunc(@(t) starfish(t,narms1,amp,[0;0.5],[],1),cparams,pref); 
% chnkr_c = chunkerfunc(@(t) starfish(t,narms2,amp,[],[],3),cparams,pref); 
% chnkr_c = reverse(chnkr_c);
chnkr = merge([chnkr_f,chnkr_c]);
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

% plot geometry and data

figure(1)
clf
plot(chnkr,'-')
hold on
quiver(chnkr)
% quiver(chnkr_f)
% quiver(chnkr_c)
axis equal

%% Get system matrix and right hand side
ibc = 1;

start = tic;
if ibc == 2
[sys,H_f] = mix_sysmat(chnkr_f,chnkr_c,zk,nu);
elseif ibc == 1
[sys,H_f] = free_sysmat(chnkr,zk,nu);
else
[sys] = clamped_sysmat(chnkr,zk);
end
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

% building RHS
src =[]; src.r = [1;2];
% src =[]; src.r = [0.5;0];src.r = [0;4.5];
free_bc_kern = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_bcs', nu);
clamp_bc_kern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_bcs', nu);

flex_free_kern = @(s,t) chnk.flex2d.kern(zk, s, t, 's');



% get free plate BCs
rhs_f = -free_bc_kern(src,chnkr_f);
% get clamped plate BCs
rhs_c = -clamp_bc_kern(src,chnkr_c);
rhs = [rhs_f;rhs_c];

if ibc == 1
rhs = -free_bc_kern(src,chnkr);
elseif ibc == 0
rhs = -clamp_bc_kern(src,chnkr);
end

%% Solve and plot
% Solving linear system

start = tic; sol = gmres(sys,rhs,[],1e-12,500); t1 = toc(start);
fprintf('%5.2e s : time for dense gmres\n',t1)    
% sol = sys\rhs;

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

start = tic; 
in = chunkerinterior(chnkr_c,{xtarg,ytarg}); out_c = ~in;
in = chunkerinterior(chnkr_f,{xtarg,ytarg}); out_f = ~in;

% in = ~chunkerinterior(reverse(chnkr_c),{xtarg,ytarg}); out_c = ~in;
out = out_c & out_f;
t1 = toc(start);

fprintf('%5.2e s : time to find points in domain\n',t1)

start1 = tic;
if ibc == 2
uscat = mixed_eval(chnkr_f,chnkr_c,sol,targets(:,out),H_f,zk,nu);
elseif ibc == 1
uscat = free_eval(chnkr,sol,targets(:,out),H_f,zk,nu);
else
uscat = clamp_eval(chnkr,sol,targets(:,out),zk);
end
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

%%
uin = flex_free_kern(src,struct("r",targets(:,out)));
utot = uscat(:)+uin(:);

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

figure(3);clf;%subplot(2,1,1)
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
title(t,"Mixed BCs")
colorbar


figure(5);clf
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
% %%
% sol_f = sol(1:2*chnkr_f.npt);
% sys_fc = free2clamp_sysmat(chnkr_f,chnkr_c,H_f,zk,nu);
% u2 = free_eval(chnkr_f,sol_f,chnkr_c,H_f,zk,nu);
% norm(u2-sys_fc(1:2:end,:)*sol_f)

function [sys,H_f] = mix_sysmat(chnkr_f,chnkr_c,zk,nu)
% build the block system matrix
    [sys_ff,H] = free_sysmat(chnkr_f,zk,nu);
    sys_fc = free2clamp_sysmat(chnkr_f,chnkr_c,H,zk,nu);
    sys_cf = clamp2free_sysmat(chnkr_f,chnkr_c,zk,nu);
    sys_cc = clamped_sysmat(chnkr_c,zk);

    sys = [sys_ff,sys_cf;sys_fc,sys_cc];
    H_f = H;
end


function [sys,H] = free_sysmat(chnkr,zk,nu)
    % assembling system matrix for free boundary

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

function sys = free2clamp_sysmat(chnkr_f,chnkr_c,H,zk,nu)
    % assemble matrix for free boundary talking to clamped boundary
    
    % defining free plate kernels
    fkern_0 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_eval', nu); 
    fkern_1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_bc1', nu); 
    
    % building system matrix
    
    sysmat0 = chunkerkernevalmat(chnkr_f,fkern_0,chnkr_c);
    sysmat1 = chunkerkernevalmat(chnkr_f,fkern_1,chnkr_c);
    
    sys = zeros(2*chnkr_c.npt,2*chnkr_f.npt);
    
    sys(1:2:end,1:2:end) = sysmat0(:,1:3:end) + sysmat0(:,2:3:end)*H;
    sys(2:2:end,1:2:end) = sysmat1(:,1:3:end) + sysmat1(:,2:3:end)*H;
    sys(1:2:end,2:2:end) = sysmat0(:,3:3:end);
    sys(2:2:end,2:2:end) = sysmat1(:,3:3:end);
end

function sys = clamp2free_sysmat(chnkr_f,chnkr_c,zk,nu)
    % assemble matrix for clamped boundary talking to free boundary

    % defining clamped plate kernels
    fkern_0 = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_bc2', nu); 
    fkern_1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_bc3', nu); 
    
    % building system matrix
    
    sysmat0 = chunkerkernevalmat(chnkr_c,fkern_0,chnkr_f);
    sysmat1 = chunkerkernevalmat(chnkr_c,fkern_1,chnkr_f);
    
    sys = zeros(2*chnkr_f.npt,2*chnkr_c.npt);
    
    sys(1:2:end,:) = sysmat0;
    sys(2:2:end,:) = sysmat1;
end

function sys = clamped_sysmat(chnkr,zk)
    % assembling system matrix for clamped boundary
    
    fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate');
    
    kappa = signed_curvature(chnkr);
    kappa = kappa(:);
    
    opts = [];
    opts.sing = 'log';
    
    sys = chunkermat(chnkr,fkern, opts);
    sys = sys - 0.5*eye(2*chnkr.npt);
    sys(2:2:end,1:2:end) = sys(2:2:end,1:2:end) + kappa.*eye(chnkr.npt);
end

function uscat = mixed_eval(chnkr_f,chnkr_c,sol,targets,H_f,zk,nu)
    % evaluate layer potentials one BC at a time
    sol_f = sol(1:2*chnkr_f.npt);
    sol_c = sol((1+2*chnkr_f.npt):end);

    uscat = free_eval(chnkr_f,sol_f,targets,H_f,zk,nu);
    uscat = uscat+clamp_eval(chnkr_c,sol_c,targets,zk);

end

function uscat = free_eval(chnkr,sol,targets,H,zk,nu)
    % evaluate layer potentials from free plate representation
    ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_eval', nu); 
    
    dens_comb = zeros(3*chnkr.npt,1);
    dens_comb(1:3:end) = sol(1:2:end);
    dens_comb(2:3:end) = H*sol(1:2:end);
    dens_comb(3:3:end) = sol(2:2:end);
    
    uscat = chunkerkerneval(chnkr, ikern,dens_comb,targets);
end

function uscat = clamp_eval(chnkr,sol,targets,zk)
    % evaluate layer potentials from clamped plate representation
    ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');   
    
    uscat = chunkerkerneval(chnkr, ikern,sol,targets);
end