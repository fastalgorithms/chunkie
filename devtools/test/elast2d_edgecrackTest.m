elast2d_edgecrackTest0();
elast2d_edgecrackTest1();
elast2d_edgecrackTest2();


function elast2d_edgecrackTest0()
%ELAST2D_EDGECRACKTEST0
%
% S-shaped edge crack in the square [-1,1]^2 with mixed boundary
% conditions: traction on the top wall, displacement (Dirichlet) on the
% bottom, right and left walls. The traction-free crack Gamma starts on
% the left wall (mouth M) and ends inside (tip T).
%
% Exact solution: Kelvinlets outside the square plus the response of a
% traction-free S-shaped crack Gamma_full in the infinite plane, where
% Gamma_full is Gamma continued smoothly to the outside of the square (it
% crosses x = -1 at M). Inside the square this is an edge-crack solution
% with an open mouth, smooth on each side of the crack at M. Its response
% is computed with the Erdogan-Gupta discretization at N = 256 (validated
% in elast2d_curvedcrackTest).
%
% Representation:
%   u = Dalt[sigma_D] + S_str[sigma_N] - D_{Gamma u Gamma_out}[phi],
% Dalt on the Dirichlet walls, the string single layer (strings along
% the outward normals) on the top wall, as in elast2d_mixedTest, and the
% crack as below. phi = u(+n side) - u(-n side), n = tangent rotated
% counterclockwise.
%
% Crack: s in [-1,1] (s = -1 at M, s = 1 at T), unknown psi = phi'(s),
% phi(s) = -int_s^1 psi (phi(T) = 0, no closure condition).
%   [-1,0]: 4 panels with 16 Gauss-Legendre nodes, the first one refined
%           dyadically L times toward M; psi at the nodes.
%   [0,1]:  tip panel, psi = rho(sig)/sqrt(1-sig) (sig in [-1,1] the
%           panel variable), rho at Nt Gauss-Jacobi(-1/2,0) nodes.
% The traction on Gamma, t(x) = -int Kt(x,gamma(s)) psi(s) ds (edge
% dislocations, Kt as in elast2d_curvedcrackTest, (x-s) Kt -> Dc I/|gamma'|),
% is collocated at the Legendre nodes and at the Gauss-Jacobi(1/2,0)
% points of the tip panel:
%   self, Legendre panel: Kt = F/(x-s) with F smooth, product integration
%       with the exact PV moments of the Legendre polynomials (2 Q_n(x)),
%       F(x,x) = Dc I/|gamma'(x)| at the coincident node;
%   self, tip panel: singularity subtraction with the closed-form
%       PV int (1-sig)^(-1/2)/(x-sig) dsig;
%   nearby panels: Gauss pieces refined geometrically toward the target;
%   other panels: the panel rules.
% The grading resolves the corner singularity at the mouth without knowing
% its exponents (Test1). Beyond about 12 levels the accuracy degrades
% (nodes within ~1e-9 of M, distances computed from absolute coordinates).
%
% The mouth: a double layer whose density does not vanish at the end of
% the arc has a point dislocation there (stress ~ 1/r). Instead the
% representation continues the crack outside the square along the tangent
% at M (Gamma_out, length 0.3) with phi_out = (phi(M) + phi'(M) l) chi(l),
% chi a smooth taper, so phi and phi' are continuous across M (phi'(M)
% from the first panel's polynomial). Gamma_out lies outside the square
% (its panels graded toward M): its only unknowns are phi(M) and phi'(M)
% and no equations are imposed on it. The left wall is split at M (a
% chunkgraph vertex, no RCIP there) and graded toward M (16 dyadic
% levels, deeper than the crack); RCIP at the four corners (nsub = 80
% for the two corners where the boundary condition changes type).
%
% NOTE (round-off): the string traction kernel vanishes identically
% between points on the same straight line, but chunkie's log-GGQ
% self-panel quadrature evaluates it at interpolated auxiliary points
% very close to the target, where (r.n)/r^2 is computed by cancellation
% from absolute coordinates (spurious entries ~ 5e-11/h, h = panel
% length). The top wall's self block is therefore set to exact zeros.
% (Dalt is integrated with the plain rule and its analytic diagonal limit
% in chunkermat, which avoids the same problem for the Dirichlet walls.)
%
% Checks (L = 8, Nt = 8..48): u at interior targets (> 0.3 from corners,
% > 0.05 from the crack) and near the mouth on both sides of the crack,
% the opening phi along Gamma (including at M), the tip coefficient
% (SIFs), and the total traction on Gamma at 41 non-collocation points.

Nts = [8 16 24 32 48];
L = 8;
prob = setupedge(true);
res = struct();
for i = 1:numel(Nts)
    r = edgesolve(prob,L,Nts(i));
    res(i).r = r;
    fprintf(['S edge crack, Kelvinlet data, L = %d, Nt = %2d (%3d crack nodes): ' ...
        'cond(A) %8.2e | err u %8.2e, near mouth %8.2e, opening %8.2e, ' ...
        'tip %8.2e, traction resid. %8.2e\n'],L,Nts(i),r.nunk,r.cond, ...
        r.erru,r.errmouth,r.errphi,r.errtip,r.res);
end

i32 = find(Nts == 32);
r = res(i32).r;
plotsol(prob,r,1,['S-shaped edge crack (red), Kelvinlet data: traction on ' ...
    'top, displacement elsewhere']);

figure(2); clf
semilogy(Nts,arrayfun(@(x) x.r.erru,res),'o-'); hold on
semilogy(Nts,arrayfun(@(x) x.r.errphi,res),'s-');
semilogy(Nts,arrayfun(@(x) x.r.errtip,res),'d-');
semilogy(Nts,arrayfun(@(x) x.r.res,res),'^-');
xlabel('N_t (tip panel nodes)'); ylabel('relative error');
legend({'u','opening','tip','traction residual'});
title(sprintf('exact solution, L = %d mouth levels',L))

assert(r.erru < 1e-11);
assert(r.errmouth < 1e-11);
assert(r.errphi < 1e-10);
assert(r.errtip < 1e-10);
assert(r.res < 1e-8);

end


function elast2d_edgecrackTest1()
%ELAST2D_EDGECRACKTEST1
%
% the same S-shaped edge crack (traction-free) with the data of the other
% (1,1)-on-top tests: traction (1,1) on the top wall, u = 0 on the bottom,
% right and left walls. No exact solution: the problem is solved on the
% mesh (L = 12 mouth levels, Nt = 32) and on the mesh with every panel of
% the walls, the crack and the extension split in two, and the two
% solutions are compared. Checks:
%
% (a) u at interior targets (> 0.3 from corners, > 0.05 from the crack),
%     u near the mouth, and the tip coefficient, mesh vs refined mesh;
%     the difference |u_c - u_f| in the domain is plotted (figure 4);
% (b) phi(M) = 0 (both crack faces are clamped at M);
% (c) the total traction on the crack at 41 non-collocation points;
% (d) the traction on the top wall, evaluated at x - d n, d = h..5h
%     (h = 1e-3), at points with |x| <= 0.7, extrapolated to d = 0.
%
% With u = 0 on the left wall the mouth cannot open, and the
% traction-free crack faces meet the clamped wall in wedges of about 41
% (above) and 139 (below) degrees: the solution has a corner singularity
% at M, resolved by the dyadic panels toward M.

L = 12; Nt = 32;
nrm = @(w) max(vecnorm(reshape(w,2,[])));
lab = {'mesh','refined mesh'};
rr = cell(2,1); probs = cell(2,1);
for i = 1:2
    probs{i} = setupedge(false,0,i);
    rr{i} = edgesolve(probs{i},L,Nt);
    fprintf(['(1,1) on top, u = 0 elsewhere, %-12s (%4d wall nodes, %3d crack nodes): ' ...
        'crack traction resid. %8.2e, |phi(M)| %8.2e (exact 0), cond %8.2e\n'], ...
        lab{i},probs{i}.cgrph.npt,rr{i}.nunk,rr{i}.res,norm(rr{i}.phiM),rr{i}.cond);
end
rc = rr{1}; r = rr{2};
du = nrm(rc.uin - r.uin)/nrm(r.uin);
dm = nrm(rc.umouth - r.umouth)/nrm(r.umouth);
dt = norm(rc.Wtip - r.Wtip)/norm(r.Wtip);
fprintf(['(1,1) on top, u = 0 elsewhere: mesh vs refined mesh: u rel. diff %8.2e, ' ...
    'near mouth %8.2e, tip %8.2e\n'],du,dm,dt);
% (d) top traction by extrapolation
xb = [linspace(-0.7,0.7,15); ones(1,15)]; nb = [0;1]*ones(1,15);
h = 1e-3; wex = [5 -10 10 -5 1];
tb = zeros(size(xb));
for id = 1:5
    t = struct('r',xb - id*h*nb,'n',nb);
    tb = tb + wex(id)*reshape(r.tfun(t),2,[]);
end
errtop = max(vecnorm(tb - [1;1]));
fprintf('(1,1) on top, refined mesh: top traction max abs. error %8.2e, tip coefficient W(1) = (%.10f, %.10f)\n', ...
    errtop,r.Wtip);
fprintf('(1,1) on top, refined mesh: u(0,0.5) = (%.12f, %.12f)\n',r.ufun([0;0.5]));

prob = probs{2};
plotsol(prob,r,3,['S-shaped edge crack (red), traction $(1,1)$ on top, ' ...
    '$u = 0$ on the other walls']);
plotselfconv(prob,rc,r,4,['$\log_{10}|u - u_{\rm ref}|/\max|u_{\rm ref}|$, ' ...
    '$u_{\rm ref}$ on the $2\times$ refined mesh']);

assert(errtop < 1e-9);
assert(du < 1e-10);
assert(dm < 1e-8);
assert(dt < 1e-10);
assert(norm(r.phiM) < 1e-6);
assert(r.res < 1e-8);

end


function elast2d_edgecrackTest2()
%ELAST2D_EDGECRACKTEST2
%
% the S-shaped edge crack with a Robin condition t + alpha*u = f
% (alpha = 1) on the bottom wall; traction on the top wall, u prescribed
% on the left and right walls. The bottom wall carries the string single
% layer (as in elast2d_robinTest): its self block is alpha*S_str (the
% string traction vanishes on a straight edge), integrated with a self
% rule split at the target (splitlogquad).
%
% (a) Kelvinlet + infinite-plane crack data (exact solution), L = 8,
%     Nt = 32: errors as in Test0, and t + alpha*u on the bottom wall,
%     evaluated at x - d n (d = h..5h, h = 1e-3) and extrapolated to d = 0;
% (b) traction (1,1) on top, f = 0 on the bottom, u = 0 on the left and
%     right walls: refinement in L (L vs L+4), phi(M) = 0, the Robin
%     condition on the bottom and the traction on the top.

alpha = 1.0;
nrm = @(w) max(vecnorm(reshape(w,2,[])));
xb = [linspace(-0.7,0.7,15); -ones(1,15)]; nb = [0;-1]*ones(1,15);
xt = [linspace(-0.7,0.7,15); ones(1,15)]; nt = [0;1]*ones(1,15);
h = 1e-3; wex = [5 -10 10 -5 1];
    function [fb,tb] = walldata(r)
        % t + alpha*u on the bottom, t on the top, extrapolated to the wall
        fb = zeros(2,15); tb = fb;
        for id = 1:5
            tq = struct('r',xb - id*h*nb,'n',nb);
            fb = fb + wex(id)*(reshape(r.tfun(tq),2,[]) + ...
                alpha*reshape(r.ufun(tq.r),2,[]));
            tq = struct('r',xt - id*h*nt,'n',nt);
            tb = tb + wex(id)*reshape(r.tfun(tq),2,[]);
        end
    end

% (a) exact solution
prob = setupedge(true,alpha);
r = edgesolve(prob,8,32);
fb = walldata(r);
fex = reshape(prob.tref(struct('r',xb,'n',nb)),2,[]) + ...
    alpha*reshape(prob.uref(xb),2,[]);
errrob = nrm(fb-fex)/nrm(fex);
fprintf(['Robin bottom, Kelvinlet data, L = 8, Nt = 32: err u %8.2e, near mouth ' ...
    '%8.2e, opening %8.2e, tip %8.2e, traction resid. %8.2e, Robin data %8.2e\n'], ...
    r.erru,r.errmouth,r.errphi,r.errtip,r.res,errrob);
assert(r.erru < 1e-11);
assert(r.errmouth < 1e-11);
assert(r.errphi < 1e-10);
assert(r.errtip < 1e-10);
assert(r.res < 1e-8);
assert(errrob < 1e-10);

% (b) (1,1) on top, homogeneous Robin bottom, clamped sides
prob = setupedge(false,alpha);
Ls = [0 4 8 12];
rr = cell(numel(Ls),1);
for i = 1:numel(Ls)
    rr{i} = edgesolve(prob,Ls(i),32);
    fprintf(['Robin bottom, (1,1) on top, L = %2d: crack traction resid. %8.2e, ' ...
        '|phi(M)| %8.2e (exact 0)\n'],Ls(i),rr{i}.res,norm(rr{i}.phiM));
end
du = zeros(1,numel(Ls)-1); dm = du; dt = du;
for i = 1:numel(Ls)-1
    r1 = rr{i}; r2 = rr{i+1};
    du(i) = nrm(r1.uin - r2.uin)/nrm(r2.uin);
    dm(i) = nrm(r1.umouth - r2.umouth)/nrm(r2.umouth);
    dt(i) = norm(r1.Wtip - r2.Wtip)/norm(r2.Wtip);
    fprintf('Robin bottom, (1,1) on top: L = %2d vs %2d: u rel. diff %8.2e, near mouth %8.2e, tip %8.2e\n', ...
        Ls(i),Ls(i+1),du(i),dm(i),dt(i));
end
r = rr{end};
[fb,tb] = walldata(r);
errbot = nrm(fb); errtop = nrm(tb - [1;1]);
fprintf(['Robin bottom, (1,1) on top, L = %d: t + alpha u on the bottom %8.2e (exact 0), ' ...
    'top traction error %8.2e, tip coefficient W(1) = (%.10f, %.10f)\n'], ...
    Ls(end),errbot,errtop,r.Wtip);
plotsol(prob,r,5,['S-shaped edge crack (red), traction $(1,1)$ on top, ' ...
    'Robin $t + u = 0$ on the bottom (green), $u = 0$ on the sides']);

assert(errbot < 1e-9);
assert(errtop < 1e-9);
assert(du(end) < 1e-7);
assert(dm(end) < 1e-5);
assert(dt(end) < 1e-8);
assert(norm(r.phiM) < 1e-5);
assert(r.res < 1e-8);

end


function prob = setupedge(exact,alpha,ref)
% geometry, kernels, square system, boundary data. exact = true: Kelvinlet
% + infinite-plane crack data; false: (1,1) on top, zero data elsewhere.
% alpha > 0: Robin condition t + alpha*u on the bottom wall (string single
% layer there), otherwise u is prescribed on the bottom wall.
% ref = 2: every panel of the walls and of the crack (edgesolve) split in
% two (the refined mesh for self-convergence checks).
if nargin < 2, alpha = 0; end
if nargin < 3, ref = 1; end
prob.ref = ref;

rng(8675309,'twister');
lam = 1.5; mu = 2.1;
nu = lam/(2*(lam+mu)); kap = 3-4*nu;
prob.lam = lam; prob.mu = mu; prob.Dc = 2*mu/(pi*(kap+1));
c = (lam+3*mu)/(2*(lam+2*mu)); prob.c = c;
prob.kd = kernel.elast2d('d',lam,mu);
prob.kdalt = kernel.elast2d('dalt',lam,mu);
prob.kdt = kernel.elast2d('dalttrac',lam,mu);
prob.ks = kernel.elast2d_string('s',lam,mu,false);
prob.kst = kernel.elast2d_string('strac',lam,mu,false);
kerns = kernel.elast2d('s',lam,mu); kernstrac = kernel.elast2d('strac',lam,mu);
Dc = prob.Dc;
kdis = kernel(); kdis.name = 'disl'; kdis.opdims = [2 2]; kdis.sing = 'smooth';
kdis.eval = @(s,t) disltrac(t.r,t.n,s.r,Dc);
prob.kdis = kdis;

% S-shaped crack: Gamma_full(s), crosses x = -1 at s = sM; Gamma = s_f in [sM,1]
sM = -0.2; a1 = 1.1/(1-sM); om = 2*pi/(1-sM);
Gf.gam = @(s) [-1 + a1*(s-sM); 0.2*sin(om*(s-sM))];
Gf.dgam = @(s) [a1 + 0*s; 0.2*om*cos(om*(s-sM))];
Gf.d2 = @(s) [0*s; -0.2*om^2*sin(om*(s-sM))];
hs = (1-sM)/2;
G.gam = @(s) Gf.gam(sM + hs*(s+1)); G.dgam = @(s) hs*Gf.dgam(sM + hs*(s+1));
G.d2 = @(s) hs^2*Gf.d2(sM + hs*(s+1));
prob.G = G; prob.Gf = Gf; prob.sM = sM; prob.hs = hs;
M = G.gam(-1); prob.M = M;
% points near the mouth, on both sides of the crack
dn = G.dgam(-1); tM = dn/norm(dn); nM = [-tM(2); tM(1)];
tm = M + [0.02 0.05 0.1].*tM + permute([1 -1],[1 3 2]).*[0.01 0.02 0.03].*nM;
prob.tmouth = reshape(tm,2,[]);

% square: edges 1 bottom, 2 right, 3 top (traction), 4 and 5 left wall
% (split at M, graded toward M)
verts = [-1 1 1 -1 M(1); -1 -1 1 1 M(2)];
prob.verts = verts;
e2v = [1 2 3 4 5; 2 3 4 5 1];
pref = []; pref.k = 16; cp = cell(5,1);
for i = 1:3, cp{i} = struct('nchmin',16*ref); end
Lgr = 16;   % chunkgraph drops splits within ~1e-6 of M: deeper breaks the 2x mesh
sp4 = [0, (1:7)/8, 1-2.^-(4:Lgr), 1]; sp5 = [0, 2.^-(Lgr:-1:4), (1:7)/8, 1];
if ref == 2, sp4 = halve(sp4); sp5 = halve(sp5); end
cp{4} = struct('tsplits',sp4(2:end-1), 'nchmin',8);
cp{5} = struct('tsplits',sp5(2:end-1), 'nchmin',8);
cgrph = chunkgraph(verts,e2v,[],cp,pref);
prob.cgrph = cgrph;
npt = cgrph.npt;
prob.bd = struct('r',cgrph.r(:,:),'n',cgrph.n(:,:));
% boundary condition per edge: 'D' displacement, 'N' traction, 'R' Robin
bc = 'DDNDD';
if alpha > 0, bc(1) = 'R'; end
prob.alpha = alpha; prob.bc = bc;
edgeid = repelem(1:5,arrayfun(@(e) e.npt,cgrph.echnks));
prob.istop = edgeid == 3;
prob.isbot = edgeid == 1 & alpha > 0;

% kernel arrays: kern(i,j) from edge j to edge i. Sources: Dalt on the
% 'D' edges, the string single layer elsewhere. Rows: (1/c) u on 'D',
% t on 'N', t + alpha*u on 'R'. Self blocks of string sources: the string
% traction vanishes identically on a straight edge, so it is dropped
% there (see the note above), leaving zero ('N') or alpha*S_str ('R').
kern(5,5) = kernel();
for i = 1:5, for j = 1:5
    if bc(j) == 'D', ku = prob.kdalt; kt = prob.kdt;
    else, ku = prob.ks; kt = mu*prob.kst; end
    if i == j && bc(j) ~= 'D', kt = kernel.zeros(2,2); end
    switch bc(i)
        case 'D', kern(i,j) = (1/c)*ku;
        case 'N', kern(i,j) = kt;
        case 'R'
            if i == j, kern(i,j) = alpha*ku; else, kern(i,j) = kt + alpha*ku; end
    end
end, end
opts = []; opts.rcip = true; opts.rcip_ignore = 5; opts.nsub_or_tol = 80;
% S_str on a straight edge: log + jump on the self panel (elast2d_robinTest)
if alpha > 0, opts.auxquads.ggqlog = splitlogquad(pref.k); end
prob.Acc = eye(2*npt) + chunkermat(cgrph,kern,opts);
kev(1,5) = kernel(); ktr(1,5) = kernel();
for j = 1:5
    if bc(j) == 'D', kev(j) = prob.kdalt; ktr(j) = prob.kdt;
    else, kev(j) = prob.ks; ktr(j) = mu*prob.kst; end
end
prob.kev = kev; prob.ktr = ktr;

% boundary data
istop2 = kron(prob.istop(:),[1;1]) == 1;
isbot2 = kron(prob.isbot(:),[1;1]) == 1;
if exact
    ns = 6; ts = 2*pi*rand(1,ns);
    src = []; src.r = 2.5*[cos(ts); sin(ts)]; src.n = randn(2,ns);
    f = randn(2*ns,1);
    prob.ukel = @(t) kerns.eval(src,struct('r',t))*f;
    prob.tkel = @(t) kernstrac.eval(src,t)*f;
    % reference: traction-free Gamma_full in the infinite plane
    cr = ibpcrack(Gf,256,Dc,true);
    cx = crk(Gf,cr.xs);
    rhof = cr.Nself\(cr.wr.*[-prob.tkel(struct('r',cx.r,'n',cx.n)); 0; 0]);
    prob.wref = @(s) reshape(cr.P(s)*rhof,2,[]);
    [chf,tf] = curvechunker(Gf,@(t) -cos(t),@(t) sin(t),@(t) cos(t),0,pi,48);
    phif = sin(tf).*prob.wref(-cos(tf));
    prob.uref = @(t) prob.ukel(t) + chunkerkerneval(chf,prob.kd,phif,t);
    % traction of the reference crack field (dislocations, Gauss-Chebyshev)
    thf = (2*(1:256)-1)*pi/(2*256); gxi = Gf.gam(cos(thf));
    prob.tref = @(t) prob.tkel(t) - (pi/256)*disltrac(t.r,t.n,gxi,Dc)*rhof;
    ud = prob.uref(prob.bd.r);
    tt = prob.tref(prob.bd);
    rhs = ud/c; rhs(istop2) = tt(istop2);
    rhs(isbot2) = tt(isbot2) + alpha*ud(isbot2);
else
    rhs = zeros(2*npt,1);
    rhs(istop2) = 1;
    prob.tkel = @(t) zeros(2*size(t.r,2),1);
end
prob.rhs = rhs;
prob.exact = exact;
end


function plotselfconv(prob,rc,rf,ifig,ttl)
% log10 of |u_c - u_f|/max|u_f| in the domain (u_c on the mesh, u_f on the
% 2x refined mesh)
nplot = 200;
[xx,yy] = meshgrid(linspace(-0.995,0.995,nplot));
tg = [xx(:).'; yy(:).'];
in = distcurve(prob.G,tg) > 0.01;
t = tg(:,in);
uf = reshape(rf.ufun(t),2,[]);
uc = reshape(rc.ufun(t),2,[]);
zz = nan(size(xx)); zz(in) = log10(max(vecnorm(uc-uf)/max(vecnorm(uf)),1e-17));
cf = prob.G.gam(linspace(-1,1,300));
figure(ifig); clf
h = pcolor(xx,yy,zz); set(h,'EdgeColor','none'); hold on
plot(prob.cgrph,'k','LineWidth',2)
plot([-1 1],[1 1],'b','LineWidth',3)
if any(prob.isbot), plot([-1 1],[-1 -1],'Color',[0 0.6 0],'LineWidth',3); end
plot(cf(1,:),cf(2,:),'r','LineWidth',2)
axis equal; xlim([-1.05 1.05]); ylim([-1.05 1.05])
set(gca,"box","off","Xtick",[],"Ytick",[]);
colorbar()
title(ttl,'Interpreter','latex','FontSize',12)
end


function out = edgesolve(prob,L,Nt)
% solve the edge-crack problem: L dyadic levels at the mouth, Nt nodes on
% the tip panel (see the description in elast2d_edgecrackTest0)
k = 16; nb = 4; sK = 0; p = 1; Lout = 0.3;
G = prob.G; Dc = prob.Dc; c = prob.c; kd = prob.kd; kdis = prob.kdis;
cgrph = prob.cgrph; npt = cgrph.npt; bd = prob.bd; M = prob.M;
k2 = @(Mx) kron(Mx,eye(2));
istop2 = kron(prob.istop(:),[1;1]) == 1;

% Legendre panels on [-1,sK], graded toward M
H = (sK+1)/nb;
br = [-1, -1 + H*2.^(-(L:-1:0)), -1 + H*(2:nb)];
if prob.ref == 2
    % every panel split in two; the tip panel [0,1] becomes a Legendre
    % panel [0,1/2] and the tip panel [1/2,1]
    br = [halve(br), (1+sK)/2]; sK = (1+sK)/2;
end
np = numel(br)-1;
pa = br(1:end-1); pb = br(2:end); hh = (pb-pa)/2; mm = (pa+pb)/2;
[xg,wg] = lege.exps(k); xg = xg(:); wg = wg(:);
U = inv(legmat(xg,k));                               % values -> coefficients
sn = mm + hh.*xg; wn = hh.*wg;                       % k x np
% tip panel [sK,1]: psi = rho(sig)/sqrt(1-sig), s = ct + ht*sig
[yt,vt] = gaussjacobi(Nt,-0.5,0);
ht = (1-sK)/2; ct = (1+sK)/2; st = ct + ht*yt; wt = ht*vt;
Vt = inv(legmat(yt,Nt)); lrowT = @(x) legmat(x,Nt)*Vt;
ss = [sn(:); st]; ww = [wn(:); wt]; Nc = numel(ss);
gs = G.gam(ss.');
[ygj,vgj] = gaussjacobi(32,-0.5,0);
[xq16,wq16] = lege.exps(16);

    function A = intK(xq)
        % int Kt(x,gamma(s)) psi(s) ds at crack parameters xq (2 rows each)
        xq = xq(:).'; cx = crk(G,xq);
        A = disltrac(cx.r,cx.n,gs,Dc).*kron(ww.',[1 1]);
        for iq = 1:numel(xq)
            xx = xq(iq); rx = cx.r(:,iq); nx = cx.n(:,iq); spx = cx.sp(iq);
            rows = 2*iq-1:2*iq;
            for j = 1:np
                cols = 2*k*(j-1)+1:2*k*j;
                if xx >= pa(j) && xx <= pb(j)
                    % self: Kt = F/(x-s), PV moments of the Legendre basis
                    tx = (xx-mm(j))/hh(j);
                    F = disltrac(rx,nx,G.gam(sn(:,j).'),Dc).* ...
                        kron((xx-sn(:,j)).',[1 1]);
                    ic = find(abs(xx-sn(:,j)) < 1e-15*max(1,abs(xx)));
                    for i = ic(:).', F(:,2*i-1:2*i) = Dc/spx*eye(2); end
                    A(rows,cols) = F.*kron((2*qleg(tx,k))*U,[1 1]);
                else
                    if xx < pa(j), e = -1; del = (pa(j)-xx)/hh(j);
                    else, e = 1; del = (xx-pb(j))/hh(j); end
                    if del < 2
                        [tp,wp] = geopieces(e,del,xq16,wq16);
                        Kp = disltrac(rx,nx,G.gam((mm(j)+hh(j)*tp).'),Dc);
                        A(rows,cols) = (Kp.*kron(hh(j)*wp.',[1 1]))* ...
                            k2(legmat(tp,k)*U);
                    end
                end
            end
            cols = 2*k*np+1:2*Nc;
            if xx >= sK
                % self, tip panel: singularity subtraction
                sx = (xx-ct)/ht; a = sqrt(1-sx);
                J = log((sqrt(2)-a)/(sqrt(2)+a))/a;
                Kt = disltrac(rx,nx,G.gam(st.'),Dc);
                corr = (Dc/spx)*(J - sum(vt./(sx-yt)));
                A(rows,cols) = ht*Kt.*kron(vt.',[1 1]) + kron(corr*lrowT(sx),eye(2));
            else
                del = (sK-xx)/ht;
                if del < 2
                    [tp,wp] = geopieces(-1,del,xq16,wq16);
                    keep = tp < 0; tp = tp(keep); wp = wp(keep).*(1-tp).^(-0.5);
                    tp = [tp; (1+ygj)/2]; wp = [wp; sqrt(1/2)*vgj];
                    Kp = disltrac(rx,nx,G.gam((ct+ht*tp).'),Dc);
                    A(rows,cols) = (Kp.*kron(ht*wp.',[1 1]))*k2(lrowT(tp));
                end
            end
        end
    end

    function P = phirows(sv)
        % phi(s) = -int_s^1 psi at crack parameters sv
        sv = sv(:); P = zeros(numel(sv),Nc);
        for iv = 1:numel(sv)
            s0 = sv(iv);
            if s0 >= sK
                P(iv,k*np+1:end) = ht*phimat((s0-ct)/ht,wl,ql,lrowT);
            else
                j = find(s0 >= pa & s0 <= pb,1);
                t0 = (s0-mm(j))/hh(j); Pl = legmat(t0,k+1);
                Rr = [1-t0, (Pl(1:k-1) - Pl(3:k+1))./(2*(1:k-1)+1)];
                P(iv,k*(j-1)+1:k*j) = -hh(j)*Rr*U;
                P(iv,k*j+1:end) = -ww(k*j+1:end).';
            end
        end
    end
[wl,ql] = gaussjacobi(Nt+2,0,0); wl = (wl.'+1)/2; ql = ql.'/2;

% Taylor data at M and the outside extension (panels graded toward M)
spM = norm(G.dgam(-1)); tM = G.dgam(-1)/spM;
Tay = zeros(p+1,Nc); Tay(1,:) = -ww.';               % phi(M) = -int psi
Tay(2,1:k) = (legmat(-1,k)*U)/spM;                   % psi(-1)/|gamma'|
Gout.gam = @(s) M - Lout*tM*(1-s); Gout.dgam = @(s) Lout*tM + 0*s; Gout.d2 = @(s) 0*[s;s];
[chout,tau] = chunkerbr(@(t) cfun(Gout,@(t) t,@(t) 1+0*t,@(t) 0*t,t), ...
    refbr([0, (1:7)/8, 1-2.^-(4:L+4), 1],prob.ref));
ell = -(1-tau)*Lout;
[chi,dchi] = taper(1-tau,p);
Ell = ell(:).^(0:p)./factorial(0:p);
dEll = [zeros(numel(ell),1), ell(:).^(0:p-1)./factorial(0:p-1)];
Pout = (Ell.*chi(:))*Tay;
dPout = (dEll.*chi(:) + Ell.*(dchi(:)*(-1/Lout)))*Tay;

% chunkers of Gamma for the displacement of the crack: the Legendre
% panels (parameter s) and the tip panel (s = 1-(1-sK)(1-t)^2, phi smooth in t)
[chg1,s1] = chunkerbr(@(t) cfun(G,@(t) t,@(t) 1+0*t,@(t) 0*t,t), br);
smt = @(t) 1-(1-sK)*(1-t).^2;
[chg2,t2] = chunkerbr(@(t) cfun(G,smt,@(t) 2*(1-sK)*(1-t),@(t) -2*(1-sK)+0*t,t), ...
    linspace(0,1,8*prob.ref+1));
Phi1 = phirows(s1); Phi2 = phirows(smt(t2));

ucrmat = @(t) chunkerkernevalmat(chg1,kd,t)*k2(Phi1) + ...
    chunkerkernevalmat(chg2,kd,t)*k2(Phi2) + chunkerkernevalmat(chout,kd,t)*k2(Pout);
tcrmat = @(t) -(disltrac(t.r,t.n,gs,Dc).*kron(ww.',[1 1])) - ...
    chunkerkernevalmat(chout,kdis,t)*k2(dPout);
tsqmat = @(t) cell2mat(arrayfun(@(j) chunkerkernevalmat(cgrph.echnks(j), ...
    prob.ktr(j),t),1:5,'UniformOutput',false));
crackrows = @(xx) -intK(xx) - chunkerkernevalmat(chout,kdis,ptinfo(G,xx))*k2(dPout);

% system: square densities and crack unknowns
xsT = gaussjacobi(Nt,0.5,0);
xcol = [sn(:); ct + ht*xsT].';
Ncr = crackrows(xcol);
Tsc = tsqmat(ptinfo(G,xcol));
cxr = crk(G,xcol);
rw = kron([ones(k*np,1); sqrt(1-xsT)].*cxr.sp(:),[1;1]);
Tsc = rw.*Tsc; Ncr = rw.*Ncr;
isbot2 = kron(prob.isbot(:),[1;1]) == 1;
isD = ~prob.istop & ~prob.isbot;
Ucs = zeros(2*npt,2*Nc);
Ucs(~istop2 & ~isbot2,:) = (1/c)*ucrmat(struct('r',bd.r(:,isD)));
Ucs(istop2,:) = tcrmat(struct('r',bd.r(:,prob.istop),'n',bd.n(:,prob.istop)));
if any(prob.isbot)
    bb = struct('r',bd.r(:,prob.isbot),'n',bd.n(:,prob.isbot));
    Ucs(isbot2,:) = tcrmat(bb) + prob.alpha*ucrmat(struct('r',bb.r));
end
AA = [prob.Acc, Ucs; Tsc, Ncr];
xsol = AA\[prob.rhs; zeros(2*Nc,1)];
sig = xsol(1:2*npt); q = xsol(2*npt+1:end);

nrm = @(w) max(vecnorm(reshape(w,2,[])));
ufun = @(t) chunkerkerneval(cgrph,prob.kev,sig,t) + ucrmat(struct('r',t))*q;
tfun = @(t) chunkerkerneval(cgrph,prob.ktr,sig,t) + tcrmat(t)*q;
out.ufun = ufun; out.tfun = tfun; out.Gout = Gout; out.nunk = Nc;
out.cond = cond(AA);
out.phiM = reshape(k2(Tay(1,:))*q,2,1);
out.Wtip = -2*sqrt(ht)*reshape(k2(lrowT(1))*q(2*k*np+1:end),2,1);   % phi ~ sqrt(1-s) W
% interior targets and points near the mouth
[xx,yy] = meshgrid(linspace(-0.98,0.98,50));
tg = [xx(:).'; yy(:).'];
tg = tg(:,min(vecnorm(tg - permute(prob.verts(:,1:4),[1 3 2]),2,1),[],3) > 0.3);
tg = tg(:,distcurve(G,tg) > 0.05);
out.uin = ufun(tg);
out.umouth = ufun(prob.tmouth);
% traction residual on the crack at 41 non-collocation points
sr = linspace(-0.98,0.98,41)+0.00731;
tres = crackrows(sr)*q + tsqmat(ptinfo(G,sr))*sig;
tscale = nrm(prob.tkel(ptinfo(G,sr)));
if tscale == 0, tscale = 1; end                    % (1,1) data: scale 1
out.res = nrm(tres)/tscale;

if prob.exact
    sv = linspace(-1,1,201);
    phin = reshape(k2(phirows(sv))*q,2,[]);
    sfs = prob.sM + prob.hs*(sv+1);
    phir = sqrt(1-sfs.^2).*prob.wref(sfs);
    out.errphi = nrm(phin-phir)/nrm(phir);
    Wtr = sqrt(1-prob.sM)*prob.wref(1);
    out.errtip = norm(out.Wtip-Wtr)/norm(Wtr);
    uex = prob.uref(tg);
    out.erru = nrm(out.uin-uex)/nrm(uex);
    out.errmouth = nrm(out.umouth-prob.uref(prob.tmouth))/nrm(uex);
end
end


function plotsol(prob,r,ifig,ttl)
nplot = 150;
[xx,yy] = meshgrid(linspace(-0.995,0.995,nplot));
tg = [xx(:).'; yy(:).'];
in = distcurve(prob.G,tg) > 0.01;
t = tg(:,in);
u = reshape(r.ufun(t),2,[]);
fplot = {u(1,:), u(2,:), vecnorm(u)};
tlabs = {'$u_1$', '$u_2$', '$|u|$'};
if prob.exact
    uex = reshape(prob.uref(t),2,[]);
    fplot{4} = log10(vecnorm(u-uex)/max(vecnorm(uex)));
    tlabs{4} = '$\log_{10}$(error)';
end
cf = prob.G.gam(linspace(-1,1,300)); co = r.Gout.gam(linspace(0,1,50));
figure(ifig); clf
nsp = numel(fplot);
if nsp == 3, set(gcf,'Position',[100 100 1200 420]); end
for ip = 1:nsp
    if nsp == 4, subplot(2,2,ip); else, subplot(1,3,ip); end
    zz = nan(size(xx)); zz(in) = fplot{ip};
    h = pcolor(xx,yy,zz); set(h,'EdgeColor','none'); hold on
    plot(prob.cgrph,'k','LineWidth',2)
    plot([-1 1],[1 1],'b','LineWidth',3)
    if any(prob.isbot), plot([-1 1],[-1 -1],'Color',[0 0.6 0],'LineWidth',3); end
    plot(cf(1,:),cf(2,:),'r','LineWidth',2)
    plot(co(1,:),co(2,:),'--','Color',[0.5 0.5 0.5],'LineWidth',1)
    axis equal; xlim([-1.3 1.05]); ylim([-1.05 1.05])
    set(gca,"box","off","Xtick",[],"Ytick",[]);
    colorbar()
    title(tlabs{ip},'Interpreter','latex','FontSize',12)
end
if exist('sgtitle','file')
    sgtitle(ttl,'Interpreter','latex')
end
end


function t = ptinfo(G,s)
c = crk(G,s); t = struct('r',c.r,'n',c.n);
end
function [ch,tn] = curvechunker(G,smap,ds,d2s,ta,tb,nch)
% chunker of gamma(smap(t)) on nch uniform panels in t; returns t at nodes
fc = @(t) cfun(G,smap,ds,d2s,t);
cparams = []; cparams.ta = ta; cparams.tb = tb; cparams.ifclosed = false;
pref = []; pref.k = 16;
ch = sort(chunkerfuncuni(fc,nch,cparams,pref));
xl = lege.exps(16);
tn = ta + (tb-ta)*((0:nch-1) + (xl(:)+1)/2)/nch; tn = tn(:).';
rr = cfun(G,smap,ds,d2s,tn);
assert(max(vecnorm(rr - ch.r(:,:))) < 1e-12,'curvechunker: node parameters');
end
function [r,d,d2] = cfun(G,smap,ds,d2s,t)
t = t(:).'; s = smap(t);
r = G.gam(s); g1 = G.dgam(s);
d = g1.*ds(t);
d2 = G.d2(s).*ds(t).^2 + g1.*d2s(t);
end
function [x,w] = gaussjacobi(n,a,b)
k = (0:n-1).'; ab = a+b;
al = (b^2-a^2)./((2*k+ab).*(2*k+ab+2));
if ab == 0, al(1) = (b-a)/(ab+2); end
kk = (1:n-1).';
be = 2./(2*kk+ab).*sqrt(kk.*(kk+a).*(kk+b).*(kk+ab)./((2*kk+ab-1).*(2*kk+ab+1)));
Jm = diag(al) + diag(be,1) + diag(be,-1);
[Vv,D] = eig(Jm); [x,i] = sort(diag(D)); Vv = Vv(:,i);
mu0 = 2^(ab+1)*gamma(a+1)*gamma(b+1)/gamma(ab+2);
w = mu0*Vv(1,:).^2; w = w(:); x = x(:);
end
function P = legmat(x,n)
x = x(:); P = zeros(numel(x),n); P(:,1) = 1;
if n > 1, P(:,2) = x; end
for k = 2:n-1, P(:,k+1) = ((2*k-1)*x.*P(:,k) - (k-1)*P(:,k-1))/k; end
end
function P = phimat(s,wl,ql,lrow)
s = s(:).'; P = 0;
for j = 1:numel(wl)
    P = P + ql(j)*lrow(1-(1-s)*wl(j)^2);
end
P = -2*sqrt(1-s(:)).*P;
end
function [chi,dchi] = taper(tp,p)
chi = 1 - betainc(tp,p+1,p+1);
dchi = -tp.^p.*(1-tp).^p/beta(p+1,p+1);
end
function d = distcurve(G,tg)
cf = G.gam(linspace(-1,1,3000));
d = min(vecnorm(permute(tg,[1 3 2]) - cf),[],2); d = d(:).';
end
function cr = ibpcrack(G,N,Dc,weighted)
th = (2*(1:N)-1)*pi/(2*N); xi = cos(th);
xs = cos((1:N-1)*pi/N);
cx = crk(G,xs); cxi = crk(G,xi);
Cau = -(pi/N)*disltrac(cx.r,cx.n,cxi.r,Dc);
Nself = [Cau; (pi*Dc/sqrt(N))*kron(ones(1,N),eye(2))];
if weighted, wr = kron([sin((1:N-1)*pi/N).*cx.sp, 1].',[1;1]); else, wr = ones(2*N,1); end
Dm = (2/N)*cos((0:N-1).*th(:)).'; Dm(1,:) = Dm(1,:)/2;
Cw = -Dm(2:N,:)./(1:N-1).';
cr.Nself = wr.*Nself; cr.wr = wr; cr.xs = xs;
cr.P = @(s) kron(umat(s,N-1)*Cw,eye(2));
end
function c = crk(G,s)
s = s(:).';
c.r = G.gam(s); d = G.dgam(s); c.sp = vecnorm(d);
c.tau = d./c.sp; c.n = [-c.tau(2,:); c.tau(1,:)];
end
function M = disltrac(xt,nt,ys,Dc)
if size(nt,2) == 1, nt = nt*ones(1,size(xt,2)); end
x = xt(1,:).' - ys(1,:); y = xt(2,:).' - ys(2,:);
r4 = (x.^2+y.^2).^2;
sxx1 = -y.*(3*x.^2+y.^2)./r4; syy1 = y.*(x.^2-y.^2)./r4; sxy1 = x.*(x.^2-y.^2)./r4;
sxx2 = x.*(x.^2-y.^2)./r4; syy2 = x.*(x.^2+3*y.^2)./r4; sxy2 = y.*(x.^2-y.^2)./r4;
n1 = nt(1,:).'; n2 = nt(2,:).';
M = zeros(2*size(xt,2),2*size(ys,2));
M(1:2:end,1:2:end) = sxx1.*n1 + sxy1.*n2; M(2:2:end,1:2:end) = sxy1.*n1 + syy1.*n2;
M(1:2:end,2:2:end) = sxx2.*n1 + sxy2.*n2; M(2:2:end,2:2:end) = sxy2.*n1 + syy2.*n2;
M = Dc*M;
end
function U = umat(s,n)
s = s(:); U = zeros(numel(s),n);
U(:,1) = 1; if n > 1, U(:,2) = 2*s; end
for k = 3:n, U(:,k) = 2*s.*U(:,k-1) - U(:,k-2); end
end

function Q = qleg(x,n)
% Q_0..Q_{n-1} (Legendre functions of the second kind), |x| < 1:
% PV int_{-1}^1 P_m(t)/(x-t) dt = 2 Q_m(x)
Q = zeros(1,n); Q(1) = 0.5*log((1+x)/(1-x));
if n > 1, Q(2) = x*Q(1) - 1; end
for m = 1:n-2, Q(m+2) = ((2*m+1)*x*Q(m+1) - m*Q(m))/(m+1); end
end

function [tp,wp] = geopieces(e,del,xq,wq)
% Gauss pieces on [-1,1] refined geometrically toward the endpoint e, for
% a target at distance del (panel units) beyond that endpoint
mmax = max(1,ceil(log2(2/del)));
d = [0, 2*2.^(-(mmax:-1:1)), 2];
tp = []; wp = [];
for i = 1:numel(d)-1
    a = d(i); b = d(i+1);
    tp = [tp; (a+b)/2 + (b-a)/2*xq(:)]; wp = [wp; (b-a)/2*wq(:)]; %#ok<AGROW>
end
tp = e - e*tp;
end

function [ch,tn] = chunkerbr(fc,brk)
% open chunker with panels [brk(i),brk(i+1)] of 16 Gauss-Legendre nodes;
% fc(t) returns r, dr/dt, d2r/dt2; tn: parameters of the nodes
k = 16; pref = []; pref.k = k;
xs = lege.exps(k);
nch = numel(brk)-1;
ch = chunker(pref); ch = ch.addchunk(nch);
tn = zeros(k,nch);
for i = 1:nch
    a = brk(i); b = brk(i+1); h = (b-a)/2;
    ts = a + (b-a)*(xs(:).'+1)/2; tn(:,i) = ts;
    [r,d,d2] = fc(ts);
    ch.rstor(:,:,i) = r; ch.dstor(:,:,i) = d*h; ch.d2stor(:,:,i) = d2*h*h;
end
adj = [(0:nch-1); (2:nch+1)]; adj(1,1) = -1; adj(2,nch) = -1;
ch.adjstor(:,1:nch) = adj;
ch.nstor(:,:,1:nch) = normals(ch);
ch.wtsstor(:,1:nch) = weights(ch);
tn = tn(:).';
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

function b = halve(b)
% breakpoints with every interval split in two
b = sort([b, (b(1:end-1)+b(2:end))/2]);
end

function b = refbr(b,ref)
if ref == 2, b = halve(b); end
end
