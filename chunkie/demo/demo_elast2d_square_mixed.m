%DEMO_ELAST2D_SQUARE_MIXED
%
% Mixed problem of 2D linear elasticity on the square [-1,1]^2: a
% cantilever, clamped on the left edge (u = 0), traction (1,1) on the top
% edge, zero traction on the bottom and right edges.
%
% Representation: alternative double layer on the Dirichlet edge, string
% single layer (strings along the outward normals) on the Neumann edges,
%
%       u = Dalt[sigma_D] + S_str[sigma_N].
%
% With c = (lam+3*mu)/(2*(lam+2*mu)) and K = 'strac' string kernel
% (traction/mu), the boundary conditions read
%
%   Dirichlet edge:  c sigma_D + Dalt sigma_D + S_str sigma_N = g
%   Neumann edges:   sigma_N + mu*K sigma_N + T[Dalt] sigma_D = t
%
% (T[Dalt] = 'dalttrac'). The Dirichlet rows are divided by c, giving the
% I + K form needed by RCIP. The operator is set up as a 4 x 4 array of
% kernels, kern(i,j) acting from edge j to edge i.
%
% At the two left corners the boundary condition changes type and the
% stress is singular; RCIP needs more levels than the default there
% (nsub = 80). Accuracy (see devtools/test/elast2d_mixedTest.m): ~6e-11
% away from the corners for a Kelvinlet exact solution.

% BEGIN DEMO ELAST2D SQUARE MIXED
%% parameters

lam = 1.5;
mu = 2.1;
c = (lam+3*mu)/(2*(lam+2*mu));

%% geometry: square as a chunkgraph (counterclockwise, outward normals)
% edges: 1 bottom, 2 right, 3 top, 4 left

verts = [-1 1 1 -1; -1 -1 1 1];
edgesendverts = [1:4; [2:4 1]];
cparams = []; cparams.nchmin = 16;
pref = []; pref.k = 16;
cgrph = chunkgraph(verts,edgesendverts,[],cparams,pref);
npt = cgrph.npt;

%% kernels: kern(i,j) from edge j (source) to edge i (target)

iD = 4;                                   % Dirichlet (clamped) edge
kd = kernel.elast2d('dalt',lam,mu);
kdt = kernel.elast2d('dalttrac',lam,mu);
ks = kernel.elast2d_string('s',lam,mu,false);
kst = kernel.elast2d_string('strac',lam,mu,false);

kern(4,4) = kernel();
for i = 1:4
    for j = 1:4
        if i == iD && j == iD
            kern(i,j) = (1/c)*kd;         % Dalt/c on the Dirichlet edge
        elseif i == iD
            kern(i,j) = (1/c)*ks;         % S_str/c, Neumann -> Dirichlet
        elseif j == iD
            kern(i,j) = kdt;              % T[Dalt], Dirichlet -> Neumann
        else
            kern(i,j) = mu*kst;           % mu*K, Neumann -> Neumann
        end
    end
end

%% boundary data: u = 0 on the left, t = (1,1) on top, t = 0 elsewhere

rhs = cell(4,1);
for i = 1:4
    ri = cgrph.echnks(i).r(:,:);
    if i == iD
        rhs{i} = zeros(2*size(ri,2),1);   % (g/c)
    else
        rhs{i} = reshape([1;1]*(abs(ri(2,:)-1) < 1e-12),[],1);
    end
end
rhs = vertcat(rhs{:});

%% build system (with RCIP at the corners) and solve

opts = []; opts.rcip = true; opts.nsub_or_tol = 80;
start = tic; sysmat = chunkermat(cgrph,kern,opts); t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1);
sys = eye(2*npt) + sysmat;

start = tic; [sigma,~,relres,iter] = gmres(sys,rhs,[],1e-11,300);
t1 = toc(start);
fprintf('%5.2e s : time for gmres, %d iters, relres %5.2e\n',t1, ...
    iter(end),relres);

%% evaluate on a grid inside the square

kev = [ks ks ks kd];                      % u = sum over edges

u0 = chunkerkerneval(cgrph,kev,sigma,[0;0]);
fprintf('u(0,0) = (%.6f, %.6f)\n',u0);

nplot = 200;
xs = linspace(-1,1,nplot+2); xs = xs(2:end-1);
[xx,yy] = meshgrid(xs,xs);
targs = [xx(:).'; yy(:).'];
in = chunkerinterior(cgrph,targs);
targin = targs(:,in);
nt = size(targin,2);

start = tic; u = chunkerkerneval(cgrph,kev,sigma,targin);
t1 = toc(start);
fprintf('%5.2e s : time to evaluate at %d targets\n',t1,nt);
u = reshape(u,2,nt);
fprintf('max |u| = %.4f\n',max(vecnorm(u)));

%% plot

fplot = {u(1,:), u(2,:), vecnorm(u)};
tlabs = {'$u_1$', '$u_2$', '$|u|$'};

figure(1); clf
tiledlayout(1,3,'TileSpacing','compact');
for ip = 1:3
    nexttile
    zz = nan(size(xx));
    zz(in) = fplot{ip};
    h = pcolor(xx,yy,zz);
    set(h,'EdgeColor','none')
    hold on
    plot(cgrph,'k','LineWidth',2)
    plot([-1 -1],[-1 1],'b','LineWidth',4)
    plot([-1 1],[1 1],'r','LineWidth',3)
    if ip == 3
        u1 = nan(size(xx)); u1(in) = u(1,:);
        u2 = nan(size(xx)); u2(in) = u(2,:);
        isk = 1:12:nplot;
        quiver(xx(isk,isk),yy(isk,isk),u1(isk,isk),u2(isk,isk),'w');
    end
    axis equal tight
    set(gca,"box","off","Xtick",[],"Ytick",[]);
    colorbar()
    title(tlabs{ip},'Interpreter','latex','FontSize',12)
end
sgtitle(['clamped left edge (blue), traction $(1,1)$ on top (red), ' ...
    'zero traction elsewhere'],'Interpreter','latex')
% END DEMO ELAST2D SQUARE MIXED
