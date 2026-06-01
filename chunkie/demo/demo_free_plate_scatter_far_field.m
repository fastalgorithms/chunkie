%DEMO_FREE_PLATE_SCATTER_FAR_FIELD
%
% Define an exterior scattering problem on a starfish-shaped domain and 
% and evaluate the far field pattern
%

% planewave vec

kvec = [8;0];
zk = norm(kvec);
nu = 1/3;

% discretize domain

cparams = [];
cparams.eps = 1.0e-6;
cparams.nover = 0;
cparams.maxchunklen = 4.0/zk; % setting a chunk length helps when the
                              % frequency is known
                              
pref = []; 
pref.k = 16;
narms = 3;
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

% defining free plate kernels

fkern1 =  @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate', nu);        % build the desired kernel
double = @(s,t) chnk.lap2d.kern(s,t,'d');
hilbert = @(s,t) chnk.lap2d.kern(s,t,'hilb');

opts = [];
opts.sing = 'log';

opts2 = [];
opts2.sing = 'pv';

% building system matrix

start = tic;
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
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

% building RHS

[~, ~, hess, third] = planewave(kvec, chnkr.r);

nx = chnkr.n(1,:).'; 
ny = chnkr.n(2,:).';

dx = chnkr.d(1,:).';
dy = chnkr.d(2,:).';

ds = sqrt(dx.*dx+dy.*dy);
taux = (dx./ds);                                                                       % normalization
tauy = (dy./ds);

kappa = signed_curvature(chnkr);
kappa = kappa(:)';

firstbc = -((hess(:,1).*(nx.*nx) + hess(:,2).*(2*nx.*ny) + hess(:,3).*(ny.*ny)) + ...
           nu.*(hess(:,1).*(taux.*taux) + hess(:,2).*(2*taux.*tauy) + hess(:,3).*(tauy.*tauy)));

secondbc = -((third(:,1).*(nx.*nx.*nx) + third(:,2).*(3*nx.*nx.*ny) +...
       third(:,3).*(3*nx.*ny.*ny) + third(:,4).*(ny.*ny.*ny))  + ...
        (2-nu).*(third(:,1).*(taux.*taux.*nx) + third(:,2).*(taux.*taux.*ny + 2*taux.*tauy.*nx) +...
        third(:,3).*(2*taux.*tauy.*ny+ tauy.*tauy.*nx) +...
        + third(:,4).*(tauy.*tauy.*ny)) + ...
        (1-nu).*kappa'.*(hess(:,1).*taux.*taux + hess(:,2).*(2*taux.*tauy) + hess(:,3).*tauy.*tauy+...
        -(hess(:,1).*nx.*nx + hess(:,2).*(2*nx.*ny) + hess(:,3).*ny.*ny)));

rhs = zeros(2*chnkr.npt, 1); rhs(1:2:end) = firstbc ; rhs(2:2:end) = secondbc;

% Solving linear system

start = tic; sol = gmres(sys,rhs,[],1e-12,500); t1 = toc(start);
fprintf('%5.2e s : time for dense gmres\n',t1)    

% evaluate at targets and plot

rmin = min(chnkr); rmax = max(chnkr);
xl = rmax(1)-rmin(1);
yl = rmax(2)-rmin(2);
nplot = 200;
xtarg = linspace(rmin(1)-xl/2,rmax(1)+xl/2,nplot); 
ytarg = linspace(rmin(2)-yl/2,rmax(2)+yl/2,nplot);
[xxtarg,yytarg] = meshgrid(xtarg,ytarg);
targets = zeros(2,length(xxtarg(:)));
targets(1,:) = xxtarg(:); targets(2,:) = yytarg(:);

start = tic; in = chunkerinterior(chnkr,{xtarg,ytarg}); t1 = toc(start);
out = ~in;

fprintf('%5.2e s : time to find points in domain\n',t1)

ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_eval', nu); 

dens_comb = zeros(3*chnkr.npt,1);
dens_comb(1:3:end) = sol(1:2:end);
dens_comb(2:3:end) = H*sol(1:2:end);
dens_comb(3:3:end) = sol(2:2:end);

start1 = tic;
uscat = chunkerkerneval(chnkr, ikern,dens_comb,targets(:,out));
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

uin = planewave(kvec,targets(:,out));
utot = uscat(:)+uin(:);

maxu = max(abs(uin(:)));

figure(2)
clf

t = tiledlayout(1,3,'TileSpacing','compact');

nexttile
zztarg = nan(size(xxtarg));
zztarg(out) = uin;
h=pcolor(xxtarg,yytarg,imag(zztarg)); shading interp
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
h=pcolor(xxtarg,yytarg,imag(zztarg)); shading interp;
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
h=pcolor(xxtarg,yytarg,imag(zztarg)); shading interp;
set(h,'EdgeColor','none')
clim([-maxu,maxu])
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
colorbar
title('$u^{\textrm{tot}}$','Interpreter','latex','FontSize',12)
title(t,"Free Plate BCs")



% evaluate the far field representation

ts = linspace(-pi,pi,300);
ts = ts(1:end-1);

targets_ff = [cos(ts); sin(ts)];

ikern_ff = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_eval_ff', nu); 

start1 = tic;
uscat = chunkerkerneval(chnkr, ikern_ff,dens_comb,targets_ff);
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

figure(3)
clf 

t = tiledlayout(1,3,'TileSpacing','compact');
title(t,'Free Plate Far Field Pattern')

nexttile
plot(ts,real(uscat));
title('Real')

nexttile
plot(ts,imag(uscat));
title('Imag')

nexttile
plot(ts,abs(uscat));
title('Abs')