elast2d_mixedTest0();
elast2d_mixedTest1();


function elast2d_mixedTest0()
%ELAST2D_MIXEDTEST
%
% analytic solution test for a mixed problem of 2D linear elasticity on
% the square [-1,1]^2: displacement (Dirichlet) data on the left edge,
% traction (Neumann) data on the other three edges. The data comes from
% a sum of Kelvinlets (point forces) placed outside the domain.
%
% Representation: alternative double layer on the Dirichlet edge, string
% single layer on the Neumann edges,
%
%       u = Dalt[sigma_D] + S_str[sigma_N],
%
% (strings along the outward normals). With c = (lam+3*mu)/(2*(lam+2*mu))
% and K = 'strac' string kernel (traction/mu), the boundary conditions
% read
%
%   Dirichlet edge:  c sigma_D + Dalt sigma_D + S_str sigma_N = g
%   Neumann edges:   sigma_N + mu*K sigma_N + T[Dalt] sigma_D = t
%
% (T[Dalt] = 'dalttrac', the traction of Dalt). The Dirichlet rows are
% divided by c so the system has the form I + K needed by RCIP. Dalt and
% K are zero on straight edges, and S_str and T[Dalt] only act between
% different edges, so the only singular interactions are at the corners,
% where RCIP handles them.
%
% Notes:
% - at the two corners where the boundary condition changes type the
%   solution has a stress singularity, and RCIP needs more levels than
%   the default to resolve it (nsub = 80).
% - with RCIP, the density on the corner panels is the compressed one,
%   which is only accurate for targets well separated from the corner on
%   the scale of the coarse panels; the error is measured away from the
%   corners (the plot shows it everywhere).

iseed = 8675309;
rng(iseed,'twister');

lam = 1.5;
mu = 2.1;
nchmin = 16;
nsub = 80;

verts = [-1 1 1 -1; -1 -1 1 1];

% point forces outside the domain

ns = 6;
ts = 2*pi*rand(1,ns);
src = []; src.r = 2.5*[cos(ts); sin(ts)]; src.n = randn(2,ns);
f = randn(2*ns,1);

kerns = kernel.elast2d('s',lam,mu);
kernstrac = kernel.elast2d('strac',lam,mu);

% boundary data: Kelvinlet displacement on the Dirichlet edge, traction
% on the Neumann edges

uex_fun = @(x) kerns.eval(src,x)*f;
tex_fun = @(x) kernstrac.eval(src,x)*f;

start = tic;
[cgrph,sigma,kev] = mixedsolve(lam,mu,nchmin,nsub,uex_fun,tex_fun);
t1 = toc(start);
fprintf('%5.2e s : time to build and solve, npt = %d\n',t1,cgrph.npt);

% targets on a grid inside the geometry

nplot = 100;
xs = linspace(-1,1,nplot+2); xs = xs(2:end-1);
[xx,yy] = meshgrid(xs,xs);
targs = [xx(:).'; yy(:).'];
in = chunkerinterior(cgrph,targs);
targin = targs(:,in);
nt = size(targin,2);

targinfo = []; targinfo.r = targin;
uex = kerns.eval(src,targinfo)*f;
usol = chunkerkerneval(cgrph,kev,sigma,targin);
usol = usol(:);

% errors

err = reshape(uex-usol,2,nt);
unorm = max(vecnorm(reshape(uex,2,nt)));
errpt = vecnorm(err)/unorm;
relerr = norm(err(:))/norm(uex);
fprintf('all targets: relative l2 error %5.2e, max pointwise rel. error %5.2e\n', ...
    relerr,max(errpt));

dcorner = min(vecnorm(targin - permute(verts,[1 3 2]),2,1),[],3);
far = dcorner > 0.3;
errfar = max(errpt(far));
fprintf('targets > 0.3 from corners: max pointwise rel. error %5.2e\n', ...
    errfar);

usol2 = reshape(usol,2,nt);
plotsol(cgrph,xx,yy,in,usol2,errpt, ...
    'mixed problem: displacement on the left, traction elsewhere');

assert(errfar < 1e-9);

end


function elast2d_mixedTest1()
%ELAST2D_MIXEDTEST1
%
% cantilever: clamped left edge (u = 0), traction (1,1) on the top edge,
% zero traction on the bottom and right edges. There is no exact
% solution, so two checks are made, away from the corners:
%
% (a) refinement: the solutions with nchmin and 2*nchmin agree on a grid
%     of interior targets more than 0.3 from the corners;
% (b) boundary condition (on the nchmin discretization): u on the left
%     edge and the traction on the other edges are evaluated at x - d*n,
%     d = h, 2h, ..., 5h (h = 1e-3), for points x on each edge with edge
%     parameter in [-0.7,0.7], extrapolated to d = 0 (error O(h^5)), and
%     compared with the data.

lam = 1.5;
mu = 2.1;
nchmin = 16;
nsub = 80;

verts = [-1 1 1 -1; -1 -1 1 1];

g_fun = @(x) zeros(2*size(x.r(:,:),2),1);
t_fun = @(x) reshape([1;1]*(abs(x.r(2,:)-1) < 1e-12),[],1);

% interior targets more than 0.3 from the corners

xs = linspace(-1,1,52); xs = xs(2:end-1);
[xx,yy] = meshgrid(xs,xs);
targs = [xx(:).'; yy(:).'];
dcorner = min(vecnorm(targs - permute(verts,[1 3 2]),2,1),[],3);
targs = targs(:,dcorner > 0.3);

nchs = [nchmin 2*nchmin];
u = cell(2,1);
for ich = 1:2
    [cgrph,sigma,kev,ktr] = mixedsolve(lam,mu,nchs(ich),nsub,g_fun,t_fun);
    u{ich} = reshape(chunkerkerneval(cgrph,kev,sigma,targs),2,[]);
    if ich == 1
        cg1 = cgrph; sig1 = sigma;
    end
end

% (a) refinement

errref = max(vecnorm(u{1}-u{2}))/max(vecnorm(u{2}));
fprintf('cantilever: nchmin %d vs %d, max rel. difference %5.2e\n', ...
    nchs(1),nchs(2),errref);

% (b) boundary condition. Edges (counterclockwise): 1 bottom, 2 right,
% 3 top, 4 left (Dirichlet).

sb = linspace(-0.7,0.7,15); e1 = ones(size(sb)); e0 = zeros(size(sb));
xb = {[sb; -e1], [e1; sb], [sb; e1], [-e1; sb]};
nb = {[e0; -e1], [e1; e0], [e0; e1], [-e1; e0]};
h = 1e-3;
wex = [5 -10 10 -5 1];
errbc = 0;
for iedge = 1:4
    bc = zeros(size(xb{iedge}));
    for id = 1:5
        t = []; t.r = xb{iedge} - id*h*nb{iedge}; t.n = nb{iedge};
        if iedge == 4
            v = chunkerkerneval(cg1,kev,sig1,t);
        else
            v = chunkerkerneval(cg1,ktr,sig1,t);
        end
        bc = bc + wex(id)*reshape(v,2,[]);
    end
    x = []; x.r = xb{iedge};
    if iedge == 4
        fb = reshape(g_fun(x),2,[]);
    else
        fb = reshape(t_fun(x),2,[]);
    end
    errbc = max(errbc,max(vecnorm(bc - fb)));
end
fprintf('cantilever: boundary condition max abs. error %5.2e\n',errbc);

% plot the solution

nplot = 100;
xs = linspace(-1,1,nplot+2); xs = xs(2:end-1);
[xx,yy] = meshgrid(xs,xs);
tp = [xx(:).'; yy(:).'];
in = chunkerinterior(cgrph,tp);
up = reshape(chunkerkerneval(cgrph,kev,sigma,tp(:,in)),2,[]);
fprintf('cantilever: u(0,0) = (%.6f, %.6f), max |u| = %.4f\n', ...
    chunkerkerneval(cgrph,kev,sigma,[0;0]),max(vecnorm(up)));

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
    plot([-1 -1],[-1 1],'b','LineWidth',4)
    plot([-1 1],[1 1],'r','LineWidth',3)
    axis equal tight
    set(gca,"box","off","Xtick",[],"Ytick",[]);
    colorbar()
    title(tlabs{ip},'Interpreter','latex','FontSize',12)
end
if exist('sgtitle','file')
    sgtitle(['clamped left (blue), traction $(1,1)$ on top (red), ' ...
        '$0$ elsewhere'],'Interpreter','latex')
end

assert(errref < 2e-9);
assert(errbc < 1e-8);

end


function [cgrph,sigma,kev,ktr] = mixedsolve(lam,mu,nchmin,nsub,g_fun,t_fun)
% solve the mixed problem on the square: displacement g on the left edge
% (edge 4), traction t on edges 1-3. g_fun and t_fun take a struct with
% fields r, n and return a 2*npt column vector. Returns the density, the
% 1 x 4 kernel array kev for evaluating u, and ktr for evaluating the
% traction (needs target normals).

c = (lam+3*mu)/(2*(lam+2*mu));
verts = [-1 1 1 -1; -1 -1 1 1];
edgesendverts = [1:4; [2:4 1]];
cparams = []; cparams.nchmin = nchmin;
pref = []; pref.k = 16;
cgrph = chunkgraph(verts,edgesendverts,[],cparams,pref);

iD = 4;
kd = kernel.elast2d('dalt',lam,mu);
kdt = kernel.elast2d('dalttrac',lam,mu);
ks = kernel.elast2d_string('s',lam,mu,false);
kst = kernel.elast2d_string('strac',lam,mu,false);

% kern(i,j): target edge i, source edge j
kern(4,4) = kernel();
for i = 1:4
    for j = 1:4
        if i == iD && j == iD
            kern(i,j) = (1/c)*kd;
        elseif i == iD
            kern(i,j) = (1/c)*ks;
        elseif j == iD
            kern(i,j) = kdt;
        else
            kern(i,j) = mu*kst;
        end
    end
end

rhs = cell(4,1);
for i = 1:4
    b = []; b.r = cgrph.echnks(i).r(:,:); b.n = cgrph.echnks(i).n(:,:);
    if i == iD
        rhs{i} = g_fun(b)/c;
    else
        rhs{i} = t_fun(b);
    end
end
rhs = vertcat(rhs{:});

opts = []; opts.rcip = true; opts.nsub_or_tol = nsub;
sys = eye(2*cgrph.npt) + chunkermat(cgrph,kern,opts);
[sigma,~,relres,iter] = gmres(sys,rhs,[],1e-11,300);
fprintf('nchmin = %2d: gmres %d iters, relres %5.2e\n',nchmin, ...
    iter(end),relres);

kev = [ks ks ks kd];
ktr = [mu*kst mu*kst mu*kst kdt];
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
    plot([-1 -1],[-1 1],'b','LineWidth',4)
    axis equal tight
    set(gca,"box","off","Xtick",[],"Ytick",[]);
    colorbar()
    title(tlabs{ip},'Interpreter','latex','FontSize',12)
end
if exist('sgtitle','file')
    sgtitle(ttl,'Interpreter','latex')
end
end
