%DEMO_CLAMPED_PLATE_SCATTER
%
% Define an exterior scattering problem on a starfish-shaped domain and 
% solve
%

% planewave vec

kvec = 10*[1.5;-1.5];
zk = norm(kvec);

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
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

% plot geometry and data

figure(1)
clf
plot(chnkr,'-x')
hold on
quiver(chnkr)
axis equal

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

% building RHS

[r1, grad] = planewave(kvec, chnkr.r);

nx = chnkr.n(1,:); 
ny = chnkr.n(2,:);

normalderiv = grad(:, 1).*(nx.')+ grad(:, 2).*(ny.');                                % Dirichlet and Neumann BC(Clamped BC)                         

firstbc = -r1;
secondbc = -normalderiv;

rhs = zeros(2*chnkr.npt, 1); rhs(1:2:end) = firstbc ; rhs(2:2:end) = secondbc;

% Solving linear system

start = tic; sol = gmres(sys,rhs,[],1e-12,100); t1 = toc(start);
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

ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');                              % build the kernel of evaluation          

start1 = tic;
uscat = chunkerkerneval(chnkr, ikern,sol, targets(:,out));
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

uin = planewave(kvec,targets(:,out));
utot = uscat(:)+uin(:);

maxu = max(abs(uin(:)));

figure(3)
clf

t = tiledlayout(1,3,'TileSpacing','compact');

nexttile
zztarg = nan(size(xxtarg));
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
zztarg = nan(size(xxtarg));
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

nexttile
zztarg = nan(size(xxtarg));
zztarg(out) = utot;
h=pcolor(xxtarg,yytarg,imag(zztarg)); h.FaceColor="interp";
set(h,'EdgeColor','none')
clim([-maxu,maxu])
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
colorbar
title('$u^{\textrm{tot}}$','Interpreter','latex','FontSize',12)
title(t,"Clamped Plate BCs")