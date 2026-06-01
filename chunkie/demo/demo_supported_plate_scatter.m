%DEMO_SUPPORTED_PLATE_SCATTER
%
% Define an exterior scattering problem on a starfish-shaped domain and 
% solve
%
% this demonstration requires advanced usage in which additional geometric
% info is passed to the kernels using chunker data. this demo also avoids
% precision loss in the kernel evaluators by using the "native" quadrature
% on certain components of the matrix

% planewave vec

kvec = 10*[1.5;-1.5];
zk = norm(kvec);
nu = 0.3;

% discretize domain

cparams = [];
cparams.eps = 1.0e-6;
cparams.nover = 0;
cparams.maxchunklen = 4.0/zk; % setting a chunk length helps when the
                              % frequency is known
                              
pref = []; 
pref.k = 16;
narms = 5;
amp = 0.25;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
chnkr = makedatarows(chnkr,2);
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

% plot geometry and data

figure(1)
clf
plot(chnkr,'-x')
hold on
quiver(chnkr)
axis equal

% calculating curvature info

kappa = signed_curvature(chnkr);
kp = arclengthder(chnkr,kappa);
kpp = arclengthder(chnkr,kp);

% supported plate kernels expect (d/ds) kappa in the first data row
% and (d^2/ds^2) kappa in the second data row

chnkr.data(1,:,:) = kp;
chnkr.data(2,:,:) = kpp;

% defining supported plate kernels

fkern1 =  @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_log',nu);           % build the desired kernel
fkern2 =  @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_smooth',nu);           % build the desired kernel

opts = [];
opts.sing = 'log';

opts2 = [];
opts2.quad = 'native';
opts2.sing = 'smooth';

% building system matrix

start = tic;
M = chunkermat(chnkr,fkern1, opts);
M2 = chunkermat(chnkr,fkern2, opts2);

c0 = (nu - 1)*(nu + 3)*(2*nu - 1)/(2*(3 - nu));

M(2:2:end,1:2:end) = M(2:2:end,1:2:end) + M2 + c0.*kappa(:).^2.*eye(chnkr.npt);
M = M - 0.5*eye(2*chnkr.npt);

sys =  M;
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

% building RHS

[val, ~, hess] = planewave(kvec, chnkr.r);

nx = chnkr.n(1,:).'; 
ny = chnkr.n(2,:).';

dx = chnkr.d(1,:).';
dy = chnkr.d(2,:).';

ds = sqrt(dx.*dx+dy.*dy);
taux = (dx./ds);                                                                       % normalization
tauy = (dy./ds);

firstbc = -val ;

secondbc = -((hess(:, 1).*(nx.*nx) + hess(:, 2).*(2*nx.*ny) + hess(:, 3).*(ny.*ny))+...
           nu*(hess(:, 1).*(taux.*taux) + hess(:, 2).*(2*taux.*tauy) + hess(:, 3).*(tauy.*tauy)));
    
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

ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_eval', nu); 

start1 = tic;
uscat = chunkerkerneval(chnkr,ikern,sol,targets(:,out));
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
h=pcolor(xxtarg,yytarg,imag(zztarg),"FaceColor","interp");
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
h=pcolor(xxtarg,yytarg,imag(zztarg),"FaceColor","interp");
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
h=pcolor(xxtarg,yytarg,imag(zztarg),"FaceColor","interp");
set(h,'EdgeColor','none')
clim([-maxu,maxu])
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
colorbar
title('$u^{\textrm{tot}}$','Interpreter','latex','FontSize',12)
title(t,"Supported Plate BCs")
