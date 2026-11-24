elast2d_curvedcrackTest0();
elast2d_curvedcrackTest1();


function elast2d_curvedcrackTest0()
%ELAST2D_CURVEDCRACKTEST0
%
% integration by parts (Erdogan-Gupta) for traction-free curved cracks in
% 2D linear elasticity, infinite plane.
%
% Crack Gamma = gamma([-1,1]), tangent tau = gamma'/|gamma'|, normal n =
% tau rotated counterclockwise. Crack field u_c = -D[phi] (D = chunkie
% 'd', the true double layer, source normal n), phi = u(+n side) - u(-n
% side) the opening. Its stress is that of a distribution of edge
% dislocations with Burgers vector b dl = -phi'(s) ds along Gamma, on any
% curve (part (a) checks this). The traction on Gamma is then a Cauchy
% singular integral of psi = phi' (parameter derivative):
%
%   t(x) = -int_{-1}^{1} Kt(x, gamma(s)) psi(s) ds,
%
% Kt = traction of the stress of a unit edge dislocation, Kt ~ Dc I /
% (|gamma'(s_x)| (s_x - s)) on the diagonal, Dc = 2 mu/(pi (kappa+1)).
% Discretization: psi = rho/sqrt(1-s^2), rho at the zeros xi_k of T_N,
% collocation at the zeros x_r of U_{N-1}:
%
%   t(x_r) = -(pi/N) sum_k Kt(gamma(x_r), gamma(xi_k)) rho(xi_k),
%
% plus int psi = 0 (the crack closes at both tips). x_r and xi_k never
% coincide, and (x_r - y) Kt(x_r, y) rho(y) is smooth, so no singular
% quadrature is needed: the rule is spectrally accurate. 'IBPw' also
% weights row r by |gamma'(x_r)| sqrt(1-x_r^2).
%
% (a) identity: traction of -D[phi] by finite differences (8th order,
%     h = 1e-3) of the 'd' field vs the dislocation integral, at targets
%     > 0.12 from a circular-arc and an S-shaped crack.
% (b) analytic test: circular-arc crack, radius R, half-angle al, under
%     equibiaxial remote tension s0. Crack-field traction t = -s0 n, and
%     (Sih, Paris & Erdogan 1962)
%       K_I  = s0 sqrt(pi R sin al) cos(al/2)/(1 + sin(al/2)^2),
%       K_II = s0 sqrt(pi R sin al) sin(al/2)/(1 + sin(al/2)^2).
%     With phi = sqrt(1-s^2) w(s), the stress intensity factors at the
%     tips s = -1, 1 are K = pi^(3/2) Dc (n.w, tau.w)/sqrt(|gamma'|) (the
%     sign of K_II depends on orientation; |K_II| is compared).

lam = 1.5;
mu = 2.1;
nu = lam/(2*(lam+mu)); kap = 3-4*nu; Dc = 2*mu/(pi*(kap+1));

% (a) dislocation identity on curved cracks

rng(1);
kd = kernel.elast2d('d',lam,mu);
ncf = 10; C = randn(2,ncf).*exp(-(0:ncf-1)/3);   % w = sum_n c_n U_n
M = 600; [sm,qm] = gc2(M);
th = (2*(1:M)-1)*pi/(2*M); xi = cos(th);
Tn = cos((1:ncf).*th(:));                         % T_{n+1}(xi)
rho = -(C.*(1:ncf))*Tn.';                         % phi' = rho/sqrt(1-s^2)
h = 1e-3; wfd = [1/280 -4/105 1/5 -4/5 0 4/5 -1/5 4/105 -1/280];
erridmax = 0;
for G = {crackgeom('arc',0.6,pi/3), crackgeom('S')}
    G = G{1}; %#ok<FXSET>
    cm = crk(G,sm); cxi = crk(G,xi);
    wm = C*umat(sm,ncf).';
    src = struct('r',cm.r,'n',cm.n);
    uc = @(t) -kd.eval(src,struct('r',t))*reshape(wm.*(qm.*cm.sp),[],1);
    tg = [0.1 0.15; 0.3 -0.2; -0.2 0.35; -0.45 -0.25; 0 0.4; 0.05 -0.1; ...
        0.6 0.1].';
    tg = tg(:,distcrack(G,tg) > 0.12);
    nt = [cos(0.3); sin(0.3)];
    err = 0; tnorm = 0;
    for it = 1:size(tg,2)
        Gu = zeros(2,2);
        for j = 1:2
            e = zeros(2,1); e(j) = 1;
            for ii = [1:4 6:9]
                Gu(:,j) = Gu(:,j) + wfd(ii)*uc(tg(:,it)+(ii-5)*h*e)/h;
            end
        end
        ep = (Gu+Gu.')/2;
        tfd = (lam*trace(ep)*eye(2) + 2*mu*ep)*nt;
        tdis = -(pi/M)*disltrac(tg(:,it),nt,cxi.r,Dc)*rho(:);
        err = max(err,norm(tfd-tdis)); tnorm = max(tnorm,norm(tfd));
    end
    fprintf('%-3s crack: dislocation identity, max rel. error %8.2e\n', ...
        G.name,err/tnorm);
    erridmax = max(erridmax,err/tnorm);
end

% (b) circular-arc crack, equibiaxial tension

R = 0.5; s0 = 1;
als = [pi/6 pi/3 2*pi/3 0.9*pi];
Ns = [8 16 32 64 128];
errK = zeros(numel(als),numel(Ns),2); condN = errK;
for ia = 1:numel(als)
    al = als(ia);
    G = crackgeom('arc',R,al);
    KIex = s0*sqrt(pi*R*sin(al))*cos(al/2)/(1+sin(al/2)^2);
    KIIex = s0*sqrt(pi*R*sin(al))*sin(al/2)/(1+sin(al/2)^2);
    ct = crk(G,[-1 1]);
    for iw = 1:2
        for iN = 1:numel(Ns)
            N = Ns(iN);
            cr = ibpcrack(G,N,Dc,iw == 2);
            cx = crk(G,cr.xs);
            rhoN = cr.Nself\(cr.wr.*[reshape(-s0*cx.n,[],1); 0; 0]);
            wt = reshape(cr.Ptip*rhoN,2,[]);
            KI = pi^1.5*Dc*sum(ct.n.*wt)./sqrt(ct.sp);
            KII = pi^1.5*Dc*sum(ct.tau.*wt)./sqrt(ct.sp);
            errK(ia,iN,iw) = max([abs(KI-KIex), abs(abs(KII)-KIIex)])/KIex;
            condN(ia,iN,iw) = cond(cr.Nself);
        end
    end
    fprintf(['arc al = %4.2f pi: K_I = %.12f, |K_II| = %.12f; rel. error ' ...
        'N=8..128 (IBPw): %s\n'],al/pi,KIex,KIIex, ...
        sprintf('%8.2e ',errK(ia,:,2)));
    fprintf(['                 cond(crack block) N=8..128: IBP %s| IBPw ' ...
        '%s\n'],sprintf('%7.1f ',condN(ia,:,1)),sprintf('%6.2f ',condN(ia,:,2)));
end

figure(3); clf
subplot(1,2,1)
for ia = 1:numel(als)
    semilogy(Ns,errK(ia,:,2),'o-'); hold on
end
xlabel('N'); ylabel('rel. error in K_I, |K_{II}|');
legend(arrayfun(@(a) sprintf('\\alpha = %.2f\\pi',a/pi),als, ...
    'UniformOutput',false))
title('arc crack, equibiaxial tension (IBPw)')
subplot(1,2,2)
for ia = 1:numel(als)
    loglog(Ns,condN(ia,:,1),'o-'); hold on
end
set(gca,'ColorOrderIndex',1)
for ia = 1:numel(als)
    loglog(Ns,condN(ia,:,2),'s--');
end
xlabel('N'); ylabel('cond(crack block)');
title('IBP (solid) vs IBPw (dashed)')

assert(erridmax < 1e-10);
i64 = find(Ns == 64);
assert(all(errK(:,i64:end,:) < 1e-11,'all'));
assert(all(abs(condN(:,end,2)./condN(:,2,2) - 1) < 1e-2));  % IBPw: flat in N
assert(all(condN(:,end,1) > 5*condN(:,2,1)));              % IBP: grows ~ N

end


function elast2d_curvedcrackTest1()
%ELAST2D_CURVEDCRACKTEST1
%
% square [-1,1]^2 with displacement (Dirichlet) data on its boundary and
% a traction-free curved crack inside: a circular arc (R = 0.6, half-angle
% pi/3) and an S-shaped crack gamma(s) = (0.5 s, 0.2 sin(pi s)).
%
% Exact solution: Kelvinlets outside the square plus the response of the
% traction-free crack to them in the infinite plane. That response is
% computed with IBPw at N = 256 (crack only; validated by the analytic
% arc test in elast2d_curvedcrackTest0), so the solves below test the
% coupling and the convergence in N.
%
% Representation u = Dalt_square[sigma] - D[phi], crack equation by
% integration by parts (see elast2d_curvedcrackTest0). The square is a
% chunkgraph with RCIP at the corners, as in elast2d_dirichletTest: the
% system is
%
%   [ I + Dalt/c (RCIP)    -D/c      ] [sigma]   [g/c]
%   [ T[Dalt] (square->Gamma)  N_IBP ] [ rho ] = [ 0 ]
%
% The RCIP corner blocks only see the crack through the right-hand-side
% like coupling term, which is smooth near the corners, and the density
% sigma is the one for coarse quadrature, so the square -> crack rows use
% the plain chunkgraph weights (the crack is > 0.4 from the boundary).
% Besides the errors in u (> 0.05 from the crack and > 0.3 from the
% corners, where the RCIP density is the compressed one), in the opening and at
% the tips, the total traction on the crack is checked at the 200 points
% cos(j pi/201), which are never collocation points (independent check of
% the boundary condition between the collocation points).

lam = 1.5;
mu = 2.1;
Ns = [8 16 32 64 128];
geoms = {crackgeom('arc',0.6,pi/3), crackgeom('S')};
mnames = {'IBP','IBPw'};

res = struct();
for ig = 1:numel(geoms)
    prob = setupsquare(lam,mu,geoms{ig});
    for iw = 1:2
        for iN = 1:numel(Ns)
            r = cracksolve(prob,iw == 2,Ns(iN));
            res(ig,iw,iN).r = r;
            fprintf(['%-3s %-4s N = %3d: cond(A) %8.2e, cond(crack) %8.2e, ' ...
                'gmres %3d | err u %8.2e, opening %8.2e, tips %8.2e, ' ...
                'traction resid. %8.2e\n'],geoms{ig}.name,mnames{iw}, ...
                Ns(iN),r.condA,r.condcr,r.iter,r.erru,r.errw,r.errtip,r.res);
        end
    end
    if ig == 2
        plotsol(prob,cracksolve(prob,true,64,true));
    end
end

figure(2); clf
mk = {'o-','s-';'o--','s--'};
subplot(1,2,1)
for ig = 1:2
    for iw = 1:2
        loglog(Ns,arrayfun(@(x) x.r.condA,squeeze(res(ig,iw,:))),mk{ig,iw}); hold on
    end
end
xlabel('N (crack nodes)'); ylabel('cond(A)');
legend({'arc IBP','arc IBPw','S IBP','S IBPw'},'Location','northwest')
title('condition number')
subplot(1,2,2)
for ig = 1:2
    for iw = 1:2
        semilogy(Ns,arrayfun(@(x) x.r.erru,squeeze(res(ig,iw,:))),mk{ig,iw}); hold on
    end
end
xlabel('N (crack nodes)'); ylabel('max rel. error in u');
legend({'arc IBP','arc IBPw','S IBP','S IBPw'})
title('accuracy')

i64 = find(Ns == 64);
for ig = 1:2
    for iw = 1:2
        r = res(ig,iw,i64).r;
        assert(r.erru < 1e-10);
        assert(r.errw < 1e-10);
        assert(r.errtip < 1e-10);
        assert(r.res < 1e-9);
    end
    assert(res(ig,2,end).r.condA < 10);                        % IBPw bounded
    assert(res(ig,1,end).r.condA > 5*res(ig,1,2).r.condA);     % IBP ~ N
end

end


function prob = setupsquare(lam,mu,G)
% square, kernels, Kelvinlets, reference crack response, boundary data

rng(8675309,'twister');
nu = lam/(2*(lam+mu)); kap = 3-4*nu;
prob.G = G; prob.lam = lam; prob.mu = mu;
prob.Dc = 2*mu/(pi*(kap+1));
prob.c = (lam+3*mu)/(2*(lam+2*mu));
prob.kd = kernel.elast2d('d',lam,mu);
prob.kdalt = kernel.elast2d('dalt',lam,mu);
prob.kdt = kernel.elast2d('dalttrac',lam,mu);
kerns = kernel.elast2d('s',lam,mu);
kernstrac = kernel.elast2d('strac',lam,mu);

% Kelvinlets outside the square
ns = 6; ts = 2*pi*rand(1,ns);
src = []; src.r = 2.5*[cos(ts); sin(ts)]; src.n = randn(2,ns);
f = randn(2*ns,1);
prob.ukel = @(t) kerns.eval(src,t)*f;
prob.tkel = @(t) kernstrac.eval(src,t)*f;

% fine crack quadrature (Gauss-Chebyshev 2nd kind, includes sqrt(1-s^2))
% for the crack field away from the crack: -D[sqrt(1-s^2) w]
M = 400; [sm,qm] = gc2(M); cm = crk(G,sm);
prob.sm = sm;
prob.csrc = struct('r',cm.r,'n',cm.n);
prob.cwts = kron(qm.*cm.sp,[1 1]);
prob.ucr = @(wm,t) -prob.kd.eval(prob.csrc,t)*(prob.cwts.'.*wm);

% reference: the traction-free crack in the infinite plane (IBPw, N = 256)
cr = ibpcrack(G,256,prob.Dc,true);
cx = crk(G,cr.xs);
rho = cr.Nself\(cr.wr.*[-prob.tkel(struct('r',cx.r,'n',cx.n)); 0; 0]);
prob.wref = cr.P(sm)*rho;
prob.wtipref = cr.Ptip*rho;

% square (counterclockwise chunkgraph), RCIP at the corners
prob.verts = [-1 1 1 -1; -1 -1 1 1];
cparams = []; cparams.nchmin = 16; pref = []; pref.k = 16;
cgrph = chunkgraph(prob.verts,[1:4; [2:4 1]],[],cparams,pref);
prob.cgrph = cgrph;
prob.bd = struct('r',cgrph.r(:,:),'n',cgrph.n(:,:));
prob.bwts = kron(cgrph.wts(:).',[1 1]);
prob.g = prob.ukel(prob.bd) + prob.ucr(prob.wref,prob.bd);
opts = []; opts.rcip = true;
prob.Acc = eye(2*cgrph.npt) + chunkermat(cgrph,(1/prob.c)*prob.kdalt,opts);
prob.uex = @(t) prob.ukel(t) + prob.ucr(prob.wref,t);
end


function r = cracksolve(prob,weighted,N,keep)
% solve the square + crack problem with N crack nodes
if nargin < 4, keep = false; end
G = prob.G; npt = prob.cgrph.npt; c = prob.c;
cr = ibpcrack(G,N,prob.Dc,weighted);
cx = crk(G,cr.xs);
wts = prob.bwts;
Tsc = cr.wr.*[prob.kdt.eval(prob.bd,struct('r',cx.r,'n',cx.n)).*wts; ...
    zeros(2,2*npt)];
Ucs = (1/c)*(-prob.kd.eval(prob.csrc,prob.bd).*prob.cwts)*cr.P(prob.sm);
A = [prob.Acc, Ucs; Tsc, cr.Nself];
rhs = [prob.g/c; zeros(2*N,1)];
x = A\rhs;
[~,~,~,it] = gmres(A,rhs,[],1e-12,500);
sig = x(1:2*npt); rho = x(2*npt+1:end);
wm = cr.P(prob.sm)*rho;

r.condA = cond(A); r.condcr = cond(cr.Nself); r.iter = it(end);
nrm = @(v) max(vecnorm(reshape(v,2,[])));
r.errw = nrm(wm - prob.wref)/nrm(prob.wref);
r.errtip = nrm(cr.Ptip*rho - prob.wtipref)/nrm(prob.wtipref);

% interior, > 0.05 from the crack and > 0.3 from the corners
[xx,yy] = meshgrid(linspace(-0.98,0.98,60));
tg = [xx(:).'; yy(:).'];
tg = tg(:,distcorner(prob.verts,tg) > 0.3);
tg = tg(:,distcrack(G,tg) > 0.05);
uex = prob.uex(struct('r',tg));
u = chunkerkerneval(prob.cgrph,prob.kdalt,sig,tg) + prob.ucr(wm,struct('r',tg));
r.erru = nrm(u-uex)/nrm(uex);

% traction residual on the crack at s_j = cos(j pi/Mr): there the
% Gauss-Chebyshev Cauchy rule on the zeros of T_Mr needs no subtraction
Mr = 201;
thr = (2*(1:Mr)-1)*pi/(2*Mr); cxr = crk(G,cos(thr));
th = (2*(1:N)-1)*pi/(2*N);
Dm = (2/N)*cos((0:N-1).*th(:)).'; Dm(1,:) = Dm(1,:)/2;
rhor = kron(cos((0:N-1).*thr(:))*Dm,eye(2))*rho;   % rho at the Mr nodes
cs = crk(G,cos((1:Mr-1)*pi/Mr));
ts = struct('r',cs.r,'n',cs.n);
tc = -(pi/Mr)*disltrac(cs.r,cs.n,cxr.r,prob.Dc)*rhor;
tsq = (prob.kdt.eval(prob.bd,ts).*wts)*sig;
r.res = nrm(tc + tsq)/nrm(prob.tkel(ts));

if keep
    r.sig = sig; r.wm = wm;
end
end


function plotsol(prob,r)
G = prob.G;
nplot = 150;
[xx,yy] = meshgrid(linspace(-0.995,0.995,nplot));
tg = [xx(:).'; yy(:).'];
in = distcrack(G,tg) > 0.02;
t = tg(:,in);
u = reshape(chunkerkerneval(prob.cgrph,prob.kdalt,r.sig,t) + ...
    prob.ucr(r.wm,struct('r',t)),2,[]);
uex = reshape(prob.uex(struct('r',t)),2,[]);
err = vecnorm(u-uex)/max(vecnorm(uex));
cf = crk(G,linspace(-1,1,400));
fplot = {u(1,:), u(2,:), vecnorm(u), log10(err)};
tlabs = {'$u_1$', '$u_2$', '$|u|$', '$\log_{10}$(error)'};
figure(1); clf
for ip = 1:4
    subplot(2,2,ip)
    zz = nan(size(xx)); zz(in) = fplot{ip};
    h = pcolor(xx,yy,zz); set(h,'EdgeColor','none'); hold on
    plot(prob.cgrph,'k','LineWidth',2)
    plot(cf.r(1,:),cf.r(2,:),'r','LineWidth',2)
    axis equal tight
    set(gca,"box","off","Xtick",[],"Ytick",[]);
    colorbar()
    title(tlabs{ip},'Interpreter','latex','FontSize',12)
end
if exist('sgtitle','file')
    sgtitle(['traction-free S-shaped crack (red) in the square, ' ...
        'IBPw, $N = 64$'],'Interpreter','latex')
end
end


function cr = ibpcrack(G,N,Dc,weighted)
% integration by parts (Erdogan-Gupta) on the crack G with N nodes.
% unknowns rho(xi_k) (interleaved 2-vectors), xi_k zeros of T_N; rows:
% traction at the zeros x_r of U_{N-1}, then int phi' = 0. P(s) maps rho
% to w(s) (phi = sqrt(1-s^2) w), Ptip = P([-1 1]).
th = (2*(1:N)-1)*pi/(2*N); xi = cos(th);
xs = cos((1:N-1)*pi/N);
cx = crk(G,xs); cxi = crk(G,xi);
Cau = -(pi/N)*disltrac(cx.r,cx.n,cxi.r,Dc);
Nself = [Cau; (pi*Dc/sqrt(N))*kron(ones(1,N),eye(2))];
if weighted
    wr = kron([sin((1:N-1)*pi/N).*cx.sp, 1].',[1;1]);
else
    wr = ones(2*N,1);
end
% rho -> T coefficients d_m; phi = -sqrt(1-s^2) sum_{m>=1} d_m U_{m-1}/m
Dm = (2/N)*cos((0:N-1).*th(:)).'; Dm(1,:) = Dm(1,:)/2;
Cw = -Dm(2:N,:)./(1:N-1).';
cr.Nself = wr.*Nself; cr.wr = wr; cr.xs = xs;
cr.P = @(s) kron(umat(s,N-1)*Cw,eye(2));
cr.Ptip = cr.P([-1 1]);
end


function G = crackgeom(name,R,al)
switch name
    case 'arc'                 % circular arc, chord centered vertically
        c0 = (R + R*cos(al))/2;
        G.gam = @(s) [R*sin(al*s); R*cos(al*s) - c0];
        G.dgam = @(s) R*al*[cos(al*s); -sin(al*s)];
    case 'S'
        G.gam = @(s) [0.5*s; 0.2*sin(pi*s)];
        G.dgam = @(s) [0.5+0*s; 0.2*pi*cos(pi*s)];
end
G.name = name;
end


function c = crk(G,s)
% points, speed, tangent and normal (tangent rotated counterclockwise)
s = s(:).';
c.r = G.gam(s); d = G.dgam(s); c.sp = vecnorm(d);
c.tau = d./c.sp; c.n = [-c.tau(2,:); c.tau(1,:)];
end


function d = distcorner(verts,tg)
d = min(vecnorm(tg - permute(verts,[1 3 2]),2,1),[],3);
end


function d = distcrack(G,tg)
cf = crk(G,linspace(-1,1,3000));
d = min(vecnorm(permute(tg,[1 3 2]) - cf.r),[],2);
d = d(:).';
end


function M = disltrac(xt,nt,ys,Dc)
% traction (normal nt at xt) of the stress of a unit edge dislocation at
% ys, t = M*b; interleaved (2 nt) x (2 ns). Plane strain, Dc = 2mu/(pi
% (kappa+1)) = mu/(2 pi (1-nu)).
if size(nt,2) == 1, nt = nt*ones(1,size(xt,2)); end
x = xt(1,:).' - ys(1,:); y = xt(2,:).' - ys(2,:);
r4 = (x.^2+y.^2).^2;
sxx1 = -y.*(3*x.^2+y.^2)./r4; syy1 = y.*(x.^2-y.^2)./r4;
sxy1 = x.*(x.^2-y.^2)./r4;                         % from b_x
sxx2 = x.*(x.^2-y.^2)./r4; syy2 = x.*(x.^2+3*y.^2)./r4;
sxy2 = y.*(x.^2-y.^2)./r4;                         % from b_y
n1 = nt(1,:).'; n2 = nt(2,:).';
M = zeros(2*size(xt,2),2*size(ys,2));
M(1:2:end,1:2:end) = sxx1.*n1 + sxy1.*n2;
M(2:2:end,1:2:end) = sxy1.*n1 + syy1.*n2;
M(1:2:end,2:2:end) = sxx2.*n1 + sxy2.*n2;
M(2:2:end,2:2:end) = sxy2.*n1 + syy2.*n2;
M = Dc*M;
end


function [s,q] = gc2(n)
% Gauss-Chebyshev 2nd kind: int sqrt(1-s^2) f(s) ds ~ sum q_j f(s_j)
th = (1:n)*pi/(n+1); s = cos(th); q = pi/(n+1)*sin(th).^2;
end


function U = umat(s,n)
% U(j,k+1) = U_k(s_j), k = 0..n-1
s = s(:); U = zeros(numel(s),n);
U(:,1) = 1; if n > 1, U(:,2) = 2*s; end
for k = 3:n, U(:,k) = 2*s.*U(:,k-1) - U(:,k-2); end
end
