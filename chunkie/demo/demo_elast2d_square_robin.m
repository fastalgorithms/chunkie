%DEMO_ELAST2D_SQUARE_ROBIN
%
% Interior Robin problem of 2D linear elasticity on the square [-1,1]^2:
%
%       t + alpha*u = f      on the boundary,
%
% where t is the traction (outward normal), with f = (1,1) on the top
% edge and f = 0 on the other three edges. For alpha > 0 the problem is
% uniquely solvable for any f (no compatibility condition, no rigid
% motions).
%
% Representation: u = S_str[sigma], the "string" single layer of
% kernel.elast2d_string, with the strings along the outward normals (the
% kernel default; never crossing the domain). Its traction has interior limit sigma + mu*K sigma
% (K = 'strac' string kernel = traction/mu, no principal value part) and
% its displacement is continuous across the boundary, so the Robin
% condition becomes
%
%       (I + mu*K + alpha*S_str) sigma = f,
%
% which is already of the form I + compact needed by RCIP.
%
% Notes:
% - the string directions must depend on the local geometry only: RCIP
%   builds the corner blocks on a copy of the geometry shifted so that the
%   corner is at the origin, so strings defined relative to a fixed point
%   (e.g. radially from the center) are wrong inside RCIP.
% - on the boundary, the atan2 term of S_str is the polar angle of
%   (x - y), which flips by pi as the target passes the source, so the
%   self-panel integrand has a jump on top of the log. The default log GGQ
%   self rules do not resolve it, so the self rule is replaced by one that
%   is split at the target (splitlogquad below).
% - S_str is not scale invariant, so RCIP needs more levels than the
%   default (nsub = 60).
% An independent Python Nystrom solver of this same equation gives
% u(0,0) = (0.2166, 0.2843), max |u| = 0.663 for the problem below
% (alpha = 1, lam = 1.5, mu = 2.1).

% BEGIN DEMO ELAST2D SQUARE ROBIN
%% parameters

lam = 1.5;
mu = 2.1;
alpha = 1.0;

%% geometry: square as a chunkgraph (counterclockwise, outward normals)

verts = [-1 1 1 -1; -1 -1 1 1];
edgesendverts = [1:4; [2:4 1]];   % edge 3 (vertex 3 -> 4) is the top
cparams = []; cparams.nchmin = 8;
pref = []; pref.k = 16;
cgrph = chunkgraph(verts,edgesendverts,[],cparams,pref);
npt = cgrph.npt;
rb = cgrph.r(:,:);
nb = cgrph.n(:,:);

%% kernels

ifimage = false;

% string single layer, strings along the outward normals (used for the
% solve and for evaluation, which must use the same strings)
kerns_str = kernel.elast2d_string('s',lam,mu,ifimage);

% Robin kernel: mu*K + alpha*S_str
krob = kernel();
krob.name = 'elasticity_string_robin';
krob.opdims = [2 2];
krob.sing = 'log';
krob.eval = @(s,t) mu*chnk.elast2d_string.kern(lam,mu,s,t,'strac', ...
    ifimage) + alpha*chnk.elast2d_string.kern(lam,mu,s,t,'s',ifimage);
krob.fmm = [];

%% build system (with RCIP at the corners)

opts = []; opts.rcip = true; opts.nsub_or_tol = 60;
opts.auxquads.ggqlog = splitlogquad(pref.k);
start = tic; sysmat = chunkermat(cgrph,krob,opts); t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1);
sys = eye(2*npt) + sysmat;

%% check against an exact solution: Kelvinlets outside the square

rng(8675309,'twister');
ns = 6;
ts = 2*pi*rand(1,ns);
src = []; src.r = 2.5*[cos(ts); sin(ts)]; src.n = randn(2,ns);
f = randn(2*ns,1);
kerns = kernel.elast2d('s',lam,mu);
kernstrac = kernel.elast2d('strac',lam,mu);

bdry = []; bdry.r = rb; bdry.n = nb; bdry.d = cgrph.d(:,:);
fex = kernstrac.eval(src,bdry)*f + alpha*kerns.eval(src,bdry)*f;
[sigex,~,relres] = gmres(sys,fex,[],1e-13,200);

ntest = 200;
tt = 2*pi*rand(1,ntest); rr = 0.95*sqrt(rand(1,ntest));
targtest = [rr.*cos(tt); rr.*sin(tt)];
tinfo = []; tinfo.r = targtest;
uex = kerns.eval(src,tinfo)*f;
utest = chunkerkerneval(cgrph,kerns_str,sigex,targtest);
fprintf('Kelvinlet check: gmres relres %5.2e, relative error %5.2e\n', ...
    relres,norm(utest(:)-uex)/norm(uex));

%% boundary data: f = (1,1) on the top edge, zero elsewhere

itop = abs(rb(2,:) - 1) < 1e-12;
frob = zeros(2,npt);
frob(:,itop) = 1;
frob = frob(:);

start = tic; [sigma,~,relres,iter] = gmres(sys,frob,[],1e-13,200);
t1 = toc(start);
fprintf('%5.2e s : time for gmres, %d iters, relres %5.2e\n',t1, ...
    iter(end),relres);

u0 = chunkerkerneval(cgrph,kerns_str,sigma,[0;0]);
fprintf('u(0,0) = (%.4f, %.4f)   [Python reference: (0.2166, 0.2843)]\n', ...
    u0);

%% evaluate on a grid inside the square

nplot = 200;
xs = linspace(-1,1,nplot+2); xs = xs(2:end-1);
[xx,yy] = meshgrid(xs,xs);
targs = [xx(:).'; yy(:).'];
in = chunkerinterior(cgrph,targs);
targin = targs(:,in);
nt = size(targin,2);

start = tic; u = chunkerkerneval(cgrph,kerns_str,sigma,targin);
t1 = toc(start);
fprintf('%5.2e s : time to evaluate at %d targets\n',t1,nt);
u = reshape(u,2,nt);
fprintf('max |u| = %.3f   [Python reference: 0.663]\n',max(vecnorm(u)));

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
sgtitle(sprintf(['Robin: $t + \\alpha u = f$, $\\alpha = %g$, ' ...
    '$f = (1,1)$ on top (red), $f = 0$ elsewhere'],alpha), ...
    'Interpreter','latex')
% END DEMO ELAST2D SQUARE ROBIN


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
