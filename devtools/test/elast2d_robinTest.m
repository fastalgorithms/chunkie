elast2d_robinTest0();
elast2d_robinTest1();


function elast2d_robinTest0()
%ELAST2D_ROBINTEST
%
% analytic solution test for the interior Robin problem of 2D linear
% elasticity on a square,
%
%       t + alpha*u = f      on the boundary,
%
% discretized as a chunkgraph with RCIP corner corrections. The boundary
% data comes from a sum of Kelvinlets (point forces) placed outside the
% domain. The solution is represented with the "string" single layer of
% kernel.elast2d_string, u = S_str[sigma], with the strings along the
% outward normals (the kernel default). Its traction has interior limit
% sigma + mu*K sigma (K = 'strac' string kernel = traction/mu, no
% principal value part) and its displacement is continuous, so we solve
%
%       (I + mu*K + alpha*S_str) sigma = f,
%
% which is already of the form I + K needed by RCIP. For alpha > 0 the
% problem is uniquely solvable, so the computed solution is compared
% directly with the exact one. The solution and the log10 error are
% plotted on a grid inside the geometry.
%
% Notes:
% - the string directions must be a function of the local geometry only:
%   RCIP builds the corner blocks on a copy of the geometry shifted so
%   that the corner is at the origin, so strings defined relative to a
%   fixed point (e.g. radially from the center) point along the edges
%   there and the solve is wrong.
% - on a straight edge the atan2 term of S_str jumps by pi as the source
%   passes the target, so the self-panel integrand is log + jump. The
%   default log GGQ self rules do not resolve the jump (error ~1e-6), so
%   the self rule is replaced by one split at the target (splitlogquad).
% - S_str is not scale invariant, so RCIP needs more levels than the
%   default (nsub = 60).
% - with RCIP, the density on the corner panels is the compressed one,
%   which is only accurate for targets well separated from the corner on
%   the scale of the coarse panels; the error is measured away from the
%   corners (the plot shows it everywhere).

iseed = 8675309;
rng(iseed,'twister');

% Lame parameters

lam = 1.5;
mu = 2.1;
alpha = 1.0;

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
src = []; src.r = 2.5*[cos(ts); sin(ts)]; src.n = randn(2,ns);
f = randn(2*ns,1);

% Kelvinlet kernels (exact solution)

kerns = kernel.elast2d('s',lam,mu);
kernstrac = kernel.elast2d('strac',lam,mu);

% string single layer (representation), strings along the outward
% normals. The same kernel must be used for the solve and for evaluation.

ifimage = false;
kerns_str = kernel.elast2d_string('s',lam,mu,ifimage);

% Robin kernel mu*K + alpha*S_str

krob = kernel();
krob.name = 'elasticity_string_robin';
krob.opdims = [2 2];
krob.sing = 'log';
krob.eval = @(s,t) mu*chnk.elast2d_string.kern(lam,mu,s,t,'strac', ...
    ifimage) + alpha*chnk.elast2d_string.kern(lam,mu,s,t,'s',ifimage);
krob.fmm = [];

% Robin data on the boundary

bdry = []; bdry.r = cgrph.r(:,:); bdry.n = cgrph.n(:,:);
bdry.d = cgrph.d(:,:);
fbdry = kernstrac.eval(src,bdry)*f + alpha*kerns.eval(src,bdry)*f;

% build system (with RCIP at the corners) and solve

opts = []; opts.rcip = true; opts.nsub_or_tol = 60;
opts.auxquads.ggqlog = splitlogquad(cgrph.echnks(1).k);
start = tic; sysmat = chunkermat(cgrph,krob,opts); t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1);

sys = eye(2*cgrph.npt) + sysmat;
start = tic; [sigma,~,relres,iter] = gmres(sys,fbdry,[],1e-13,200);
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

start = tic; usol = chunkerkerneval(cgrph,kerns_str,sigma,targin); t1 = toc(start);
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
    'interior Robin problem on a square (strings)');

assert(errfar < 1e-8);

end



function elast2d_robinTest1()
%ELAST2D_ROBINTEST1
%
% the same Robin problem with the data of the demo: f = (1,1)
% on the top edge, 0 on the other three edges. There is no exact
% solution, so two checks are made, both away from the corners (where
% the data is discontinuous and the RCIP density is the compressed one):
%
% (a) refinement: the solutions with nchmin = 8 and nchmin = 32 agree on
%     a grid of interior targets more than 0.3 from the corners;
% (b) boundary condition (nchmin = 8): t + alpha*u is evaluated at x - d*n,
%     d = h, 2h, ..., 5h (h = 1e-3), for points x on each edge with
%     edge parameter in [-0.7,0.7], extrapolated to d = 0 (error O(h^5)),
%     and compared with the data.

lam = 1.5;
mu = 2.1;
alpha = 1.0;
kerns_str = kernel.elast2d_string('s',lam,mu,false);
kernstrac_str = kernel.elast2d_string('strac',lam,mu,false);
krob = kernel();
krob.name = 'elasticity_string_robin';
krob.opdims = [2 2];
krob.sing = 'log';
krob.eval = @(s,t) mu*chnk.elast2d_string.kern(lam,mu,s,t,'strac', ...
    false) + alpha*chnk.elast2d_string.kern(lam,mu,s,t,'s',false);
krob.fmm = [];
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
    opts = []; opts.rcip = true; opts.nsub_or_tol = 60;
    opts.auxquads.ggqlog = splitlogquad(pref.k);
    sys = eye(2*cgrph.npt) + chunkermat(cgrph,krob,opts);
    [sigma,~,relres,iter] = gmres(sys,rhs(:),[],1e-13,200);
    fprintf('top (1,1) data, nchmin = %2d: gmres %d iters, relres %5.2e\n', ...
        nchs(ich),iter(end),relres);
    ui = reshape(chunkerkerneval(cgrph,kerns_str,sigma,targs),2,[]);

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
    v = mu*chunkerkerneval(cg1,kernstrac_str,sig1,t) + ...
        alpha*chunkerkerneval(cg1,kerns_str,sig1,t);
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
up = reshape(chunkerkerneval(cgrph,kerns_str,sigma,tp(:,in)),2,[]);

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
    sgtitle('Robin $f = (1,1)$ on top (red), $0$ elsewhere, $\alpha = 1$','Interpreter','latex')
end

assert(errref < 1e-7);
assert(errbc < 3e-7);

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


function auxquads = splitlogquad(k)
% chunkie's ggq log rule with the self rule replaced by one that is split
% at each target node s_j, so that
%       f0(s) + log|s - s_j| f1(s) + sign(s - s_j) f2(s)
% (f0, f1, f2 smooth) is integrated accurately on the self panel. On each
% side of s_j: composite Gauss-Legendre on intervals graded geometrically
% (ratio 4) toward s_j, and on the last interval [0,e] next to s_j
% (e < 2e-3) a rule on Gauss-Legendre nodes whose weights integrate u^i
% and u^i log(u), i < m/2, exactly. All nodes stay > 1e-6 from s_j, which
% matters for 'strac', which is 0/0 on a straight edge.
auxquads = chnk.quadggq.setup(k,'log');
m = 16; q = 4; emax = 2e-3;
[xg,wg] = lege.exps(m); xg = (xg+1)/2; wg = wg/2;
ii = 0:m/2-1;
wlog = [xg.^ii, (xg.^ii).*log(xg)].' \ [1./(ii+1), -1./(ii+1).^2].';
[xleg,~] = lege.exps(k);
for j = 1:k
    [ul,wl] = graded(xleg(j)+1);
    [ur,wr] = graded(1-xleg(j));
    xs = [xleg(j) - ul; xleg(j) + ur];
    auxquads.xs0{j} = xs;
    auxquads.wts0{j} = [wl; wr];
    auxquads.ainterps0{j} = lege.matrin(k,xs);
end

    function [u,w] = graded(h)
        nlev = max(0,ceil(log(h/emax)/log(q)));
        ends = h*q.^-(0:nlev);
        a = ends(2:end); b = ends(1:end-1);
        e = ends(end);
        u = [reshape(a + (b-a).*xg, [], 1); e*xg];
        w = [reshape((b-a).*wg, [], 1); e*wlog];
    end
end
