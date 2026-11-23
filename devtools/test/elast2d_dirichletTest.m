elast2d_dirichletTest0();
elast2d_dirichletTest1();


function elast2d_dirichletTest0()
%ELAST2D_DIRICHLETTEST
%
% analytic solution test for the interior displacement (Dirichlet)
% problem of 2D linear elasticity on a square, discretized as a
% chunkgraph with RCIP corner corrections. The boundary data is the
% displacement of a sum of Kelvinlets (point forces) placed outside the
% domain. The solution is represented with the alternative (principal
% value free) double layer, u = Dalt[sigma], whose interior limit is
%
%       c sigma + Dalt sigma,    c = (lam+3*mu)/(2*(lam+2*mu)).
%
% RCIP in chunkermat assumes a system of the form I + K, so we solve
%
%       (I + Dalt/c) sigma = g/c.
%
% The Dirichlet problem is uniquely solvable, so the computed solution is
% compared directly with the exact one. The solution and the log10 error
% are plotted on a grid inside the geometry.
%
% With RCIP, the density on the corner panels is the compressed one,
% which is only accurate for targets well separated from the corner on
% the scale of the coarse panels, so the error is measured (and the
% assert applied) on targets more than 0.3 from the corners; the plot
% shows the error everywhere.

iseed = 8675309;
rng(iseed,'twister');

% Lame parameters

lam = 1.5;
mu = 2.1;
c = (lam+3*mu)/(2*(lam+2*mu));

% geometry: square [-1,1]^2 as a chunkgraph (counterclockwise, so the
% normals point out of the domain)

verts = [-1 1 1 -1; -1 -1 1 1];
edgesendverts = [1:4; [2:4 1]];
cparams = []; cparams.nchmin = 8;
pref = []; pref.k = 16;
start = tic; cgrph = chunkgraph(verts,edgesendverts,[],cparams,pref);
t1 = toc(start);
fprintf('%5.2e s : time to build geo, nch = %d, npt = %d\n',t1, ...
    cgrph.nch,cgrph.npt);

% point forces outside the domain

ns = 6;
ts = 2*pi*rand(1,ns);
src = []; src.r = 2.5*[cos(ts); sin(ts)];
f = randn(2*ns,1);

kerns = kernel.elast2d('s',lam,mu);
kernd = kernel.elast2d('dalt',lam,mu);

% displacement data on the boundary

bdry = []; bdry.r = cgrph.r(:,:);
ubdry = kerns.eval(src,bdry)*f;

% build system (with RCIP at the corners) and solve

opts = []; opts.rcip = true;
start = tic; sysmat = chunkermat(cgrph,(1/c)*kernd,opts); t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1);

sys = eye(2*cgrph.npt) + sysmat;
start = tic; [sigma,~,relres,iter] = gmres(sys,ubdry/c,[],1e-13,200);
t1 = toc(start);
fprintf('%5.2e s : time for gmres, %d iters, relres %5.2e\n',t1, ...
    iter(end),relres);

% targets on a grid inside the geometry (boundary excluded)

nplot = 100;
xs = linspace(-1,1,nplot+2); xs = xs(2:end-1);
[xx,yy] = meshgrid(xs,xs);
targs = [xx(:).'; yy(:).'];
in = chunkerinterior(cgrph,targs);
targin = targs(:,in);
nt = size(targin,2);

% exact and computed displacements

targinfo = []; targinfo.r = targin;
uex = kerns.eval(src,targinfo)*f;

start = tic; usol = chunkerkerneval(cgrph,kernd,sigma,targin); t1 = toc(start);
fprintf('%5.2e s : time to evaluate at %d targets\n',t1,nt);
usol = usol(:);

% errors

err = reshape(uex-usol,2,nt);
unorm = max(vecnorm(reshape(uex,2,nt)));
errpt = vecnorm(err)/unorm;
relerr = norm(err(:))/norm(uex);
fprintf('all targets: relative l2 error %5.2e, max pointwise rel. error %5.2e\n', ...
    relerr,max(errpt));

% targets away from the corners
dcorner = min(vecnorm(targin - permute(verts,[1 3 2]),2,1),[],3);
far = dcorner > 0.3;
errfar = max(errpt(far));
fprintf('targets > 0.3 from corners: max pointwise rel. error %5.2e\n', ...
    errfar);

% plot the computed solution and the log10 error

usol2 = reshape(usol,2,nt);
plotsol(cgrph,xx,yy,in,usol2,errpt, ...
    'interior displacement problem on a square (Dalt)');

assert(errfar < 1e-10);

end



function elast2d_dirichletTest1()
%ELAST2D_DIRICHLETTEST1
%
% the same displacement problem with the data of the demo: displacement (1,1)
% on the top edge, 0 on the other three edges. There is no exact
% solution, so two checks are made, both away from the corners (where
% the data is discontinuous and the RCIP density is the compressed one):
%
% (a) refinement: the solutions with nchmin = 8 and nchmin = 32 agree on
%     a grid of interior targets more than 0.3 from the corners;
% (b) boundary condition (nchmin = 8): u is evaluated at x - d*n,
%     d = h, 2h, ..., 5h (h = 1e-3), for points x on each edge with
%     edge parameter in [-0.7,0.7], extrapolated to d = 0 (error O(h^5)),
%     and compared with the data.

lam = 1.5;
mu = 2.1;
c = (lam+3*mu)/(2*(lam+2*mu));
kernd = kernel.elast2d('dalt',lam,mu);
verts = [-1 1 1 -1; -1 -1 1 1];
edgesendverts = [1:4; [2:4 1]];

dat = @(x) [1;1]*(abs(x(2,:)-1) < 1e-12);

% interior targets more than 0.3 from the corners

xs = linspace(-1,1,52); xs = xs(2:end-1);
[xx,yy] = meshgrid(xs,xs);
targs = [xx(:).'; yy(:).'];
dcorner = min(vecnorm(targs - permute(verts,[1 3 2]),2,1),[],3);
targs = targs(:,dcorner > 0.3);

% boundary check points and inward offsets

sb = linspace(-0.7,0.7,15); e1 = ones(size(sb)); e0 = zeros(size(sb));
xb = [sb, e1, sb, -e1; -e1, sb, e1, sb];
nb = [e0, e1, e0, -e1; -e1, e0, e1, e0];
fb = dat(xb);
h = 1e-3;
wex = [5 -10 10 -5 1];

nchs = [8 32];
u = cell(2,1);
for ich = 1:2
    cparams = []; cparams.nchmin = nchs(ich);
    pref = []; pref.k = 16;
    cgrph = chunkgraph(verts,edgesendverts,[],cparams,pref);
    rhs = dat(cgrph.r(:,:));
    opts = []; opts.rcip = true;
    sys = eye(2*cgrph.npt) + chunkermat(cgrph,(1/c)*kernd,opts);
    [sigma,~,relres,iter] = gmres(sys,rhs(:)/c,[],1e-13,200);
    fprintf('top (1,1) data, nchmin = %2d: gmres %d iters, relres %5.2e\n', ...
        nchs(ich),iter(end),relres);
    ui = reshape(chunkerkerneval(cgrph,kernd,sigma,targs),2,[]);

    u{ich} = ui;
    if ich == 1
        cg1 = cgrph; sig1 = sigma;
    end
end

% (a) refinement

errref = max(vecnorm(u{1}-u{2}))/max(vecnorm(u{2}));
fprintf('top (1,1) data: nchmin 8 vs 32, max rel. difference %5.2e\n', ...
    errref);

% (b) boundary condition (on the nchmin = 8 discretization)

bc = zeros(size(xb));
for id = 1:5
    t = []; t.r = xb - id*h*nb; t.n = nb;
    v = chunkerkerneval(cg1,kernd,sig1,t);
    bc = bc + wex(id)*reshape(v,2,[]);
end
errbc = max(vecnorm(bc - fb))/max(vecnorm(fb));
fprintf('top (1,1) data: boundary condition max rel. error %5.2e\n', ...
    errbc);

% plot the solution

nplot = 100;
xs = linspace(-1,1,nplot+2); xs = xs(2:end-1);
[xx,yy] = meshgrid(xs,xs);
tp = [xx(:).'; yy(:).'];
in = chunkerinterior(cgrph,tp);
up = reshape(chunkerkerneval(cgrph,kernd,sigma,tp(:,in)),2,[]);

fplot = {up(1,:), up(2,:), vecnorm(up)};
tlabs = {'$u_1$', '$u_2$', '$|u|$'};
figure(2); clf
for ip = 1:3
    subplot(1,3,ip)
    zz = nan(size(xx));
    zz(in) = fplot{ip};
    hh = pcolor(xx,yy,zz);
    set(hh,'EdgeColor','none')
    hold on
    plot(cgrph,'k','LineWidth',2)
    plot([-1 1],[1 1],'r','LineWidth',3)
    axis equal tight
    set(gca,"box","off","Xtick",[],"Ytick",[]);
    colorbar()
    title(tlabs{ip},'Interpreter','latex','FontSize',12)
end
if exist('sgtitle','file')
    sgtitle('displacement $(1,1)$ on top (red), $0$ elsewhere','Interpreter','latex')
end

assert(errref < 1e-9);
assert(errbc < 5e-9);

end

function plotsol(cgrph,xx,yy,in,u,errpt,ttl)
fplot = {u(1,:), u(2,:), vecnorm(u), log10(errpt)};
tlabs = {'$u_1$', '$u_2$', '$|u|$', '$\log_{10}$(error)'};
figure(1); clf
for ip = 1:4
    subplot(2,2,ip)
    zz = nan(size(xx));
    zz(in) = fplot{ip};
    h = pcolor(xx,yy,zz);
    set(h,'EdgeColor','none')
    hold on
    plot(cgrph,'k','LineWidth',2)
    axis equal tight
    set(gca,"box","off","Xtick",[],"Ytick",[]);
    colorbar()
    title(tlabs{ip},'Interpreter','latex','FontSize',12)
end
if exist('sgtitle','file')
    sgtitle(ttl,'Interpreter','latex')
end
end
