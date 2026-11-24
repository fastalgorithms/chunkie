elast2d_crackTest0();


function elast2d_crackTest0()
%ELAST2D_CRACKTEST
%
% analytic solution test for 2D linear elasticity in the unit disk with a
% traction-free straight slit Gamma = (-a,0)-(a,0), displacement
% (Dirichlet) data on the circle, and a comparison of three
% discretizations of the crack equation.
%
% Exact solution: Kelvinlets (point forces) outside the disk plus the
% field they induce on a traction-free slit in the infinite plane. On a
% straight slit the crack equation is diagonal in Chebyshev polynomials,
% so this reference is computed spectrally (to machine precision) and
% independently of the solvers below.
%
% Representation:
%   u = Dalt_circle[sigma] + u_c[phi],   u_c[phi] = -D[phi]  (D = chunkie
%   'd', the true double layer, source normal (0,1)),
% where phi = u(above) - u(below) is the crack opening. -D[phi] has jump
% phi across Gamma and continuous traction (Somigliana), so it is a crack
% field; its stress is that of a distribution of edge dislocations with
% density b = -phi' (Burgers vector per length). Both facts were checked
% numerically against finite differences of the 'd' field. On the slit
% the traction is
%
%   t(x) = -Dc * PV int phi'(y)/(x-y) dy,  Dc = 2 mu/(pi (kappa+1)),
%
% kappa = 3 - 4 nu, and with phi = a sqrt(1-s^2) w(s), w = sum c_n U_n(s),
% t = -pi Dc sum (n+1) c_n U_n(s) (n = 0 is the Griffith crack).
%
% Equations: on the circle (I + Dalt/c) sigma + u_c[phi]/c = g/c (as in
% the Dirichlet test); on the slit t(sigma) + t(phi) = 0. Crack
% formulations compared:
%
%  'H'    hypersingular: unknown w at the zeros of U_N, collocation there.
%         Order +1 operator, cond ~ N.
%  'IBP'  integration by parts (Erdogan-Gupta): unknown rho, phi' =
%         rho/sqrt(1-s^2), at the zeros of T_N; Cauchy kernel collocated
%         at the zeros of U_{N-1}, plus int phi' = 0.
%  'IBPw' same, with the collocation rows weighted by sqrt(1-x^2) (a
%         diagonal scaling under which the discrete Cauchy operator is an
%         isometry).
%  'CP'   Calderon preconditioning of 'H' from the right: w = S_omega[mu],
%         S_omega[mu] = int log|s-y| mu(y)/sqrt(1-y^2) dy, so that
%         N_omega S_omega = pi^2 I + compact (Bruno-Lintner).
%
% For each method and N the condition numbers of the full system and of
% the crack block, the GMRES iterations (tol 1e-12), and the errors in u
% (interior, > 0.05 from the slit), in the opening and at the tips
% (proportional to the stress intensity factors) are reported.

lam = 1.5;
mu = 2.1;
a = 0.5;
methods = {'H','IBP','IBPw','CP'};
Ns = [8 16 32 64 128];

prob = setup(lam,mu,a);

res = struct();
for im = 1:numel(methods)
    for iN = 1:numel(Ns)
        r = cracksolve(prob,methods{im},Ns(iN));
        res(im,iN).r = r;
        fprintf(['%-4s N = %3d: cond(A) %8.2e, cond(crack) %8.2e, ' ...
            'gmres %3d | err u %8.2e, opening %8.2e, tips %8.2e\n'], ...
            methods{im},Ns(iN),r.condA,r.condcr,r.iter,r.erru, ...
            r.errw,r.errtip);
    end
end

% plot the solution (IBPw, N = 32) and the error

r = cracksolve(prob,'IBPw',32,true);
plotsol(r,a);

% convergence and conditioning

figure(2); clf
mk = {'o-','s-','d-','^-'};
subplot(1,2,1)
for im = 1:numel(methods)
    loglog(Ns,arrayfun(@(x) x.r.condA,res(im,:)),mk{im}); hold on
end
xlabel('N (crack nodes)'); ylabel('cond(A)'); legend(methods,'Location','northwest')
title('condition number')
subplot(1,2,2)
for im = 1:numel(methods)
    semilogy(Ns,arrayfun(@(x) x.r.erru,res(im,:)),mk{im}); hold on
end
xlabel('N (crack nodes)'); ylabel('max rel. error in u'); legend(methods)
title('accuracy')

% asserts (N = 32 and the conditioning trends)

i32 = find(Ns == 32);
for im = 1:numel(methods)
    r = res(im,i32).r;
    assert(r.erru < 1e-10);
    assert(r.errw < 1e-10);
    assert(r.errtip < 1e-10);
end
condA = @(m,iN) res(strcmp(methods,m),iN).r.condA;
assert(condA('IBPw',numel(Ns)) < 10);            % bounded
assert(condA('CP',numel(Ns)) < 50);
assert(condA('H',numel(Ns)) > 10*condA('H',1));  % grows ~ N

end


function prob = setup(lam,mu,a)
% geometry, kernels, exact solution and boundary data

rng(8675309,'twister');
nu = lam/(2*(lam+mu)); kap = 3-4*nu;
prob.lam = lam; prob.mu = mu; prob.a = a;
prob.Dc = 2*mu/(pi*(kap+1));
prob.c = (lam+3*mu)/(2*(lam+2*mu));
prob.kd = kernel.elast2d('d',lam,mu);
prob.kdalt = kernel.elast2d('dalt',lam,mu);
prob.kdt = kernel.elast2d('dalttrac',lam,mu);
kerns = kernel.elast2d('s',lam,mu);
kernstrac = kernel.elast2d('strac',lam,mu);

% Kelvinlets outside the disk
ns = 6; ts = 2*pi*rand(1,ns);
src = []; src.r = 2.5*[cos(ts); sin(ts)]; src.n = randn(2,ns);
f = randn(2*ns,1);

% reference crack opening: N w = -t_Kelvin on the slit, solved
% spectrally, c_n = t_n/(pi Dc (n+1)) with t_Kelvin = sum t_n U_n
Ne = 80;
[se,qe] = gc2(Ne);
tK = reshape(kernstrac.eval(src,crk(se,a))*f,2,[]);
tKn = (2/pi)*(tK.*qe)*umat(se,Ne);
cex = tKn./(pi*prob.Dc*(1:Ne));
prob.wex = @(s) cex*umat(s,Ne).';

% fine crack quadrature (for the crack -> circle map and evaluation)
prob.M = 400;
[prob.sm,prob.qm] = gc2(prob.M);

% circle
cparams = []; cparams.maxchunklen = 0.5; pref = []; pref.k = 16;
prob.chnkr = chunkerfunc(@(t) chnk.curves.bymode(t,1),cparams,pref);
circ = []; circ.r = prob.chnkr.r(:,:); circ.n = prob.chnkr.n(:,:);
prob.circ = circ;

% exact solution and data
prob.uex = @(t) kerns.eval(src,t)*f + ucrack(prob,prob.wex(prob.sm),t);
prob.g = prob.uex(circ);

% circle block
prob.Acc = eye(2*prob.chnkr.npt) + chunkermat(prob.chnkr, ...
    (1/prob.c)*prob.kdalt);
end


function r = cracksolve(prob,method,N,keep)
% assemble and solve the coupled system for one crack formulation

if nargin < 4, keep = false; end
a = prob.a; Dc = prob.Dc; c = prob.c; M = prob.M;
npt = prob.chnkr.npt;

switch method
    case {'H','CP'}
        [sj,~] = gc2(N); xs = sj;
        V = umat(sj,N);
        P = umat(prob.sm,N)/V;                 % w(sm) from w(sj)
        Ptip = umat([-1 1],N)/V;
        Nself = -pi*Dc*V*diag(1:N)/V;          % traction at sj from w(sj)
        extra = 0;
    case {'IBP','IBPw'}
        th = (2*(1:N)-1)*pi/(2*N); xi = cos(th);
        xs = cos((1:N-1)*pi/N);
        Cau = -Dc*(pi/N)./(xs(:) - xi(:).');
        Nself = [Cau; pi*Dc/sqrt(N)*ones(1,N)];   % last row: int phi' = 0
        Tm = cos((0:N-1).*th(:));
        Dm = (2/N)*Tm.'; Dm(1,:) = Dm(1,:)/2;
        Cw = -Dm(2:N,:)./(1:N-1).';            % rho -> U coeffs of w
        P = umat(prob.sm,N-1)*Cw;
        Ptip = umat([-1 1],N-1)*Cw;
        extra = 1;
end

% circle -> crack traction (normal (0,1)) at the collocation points,
% smooth rule (the slit is 0.5 away from the circle)
wts = prob.chnkr.wts(:);
Tsc = prob.kdt.eval(prob.circ,crk(xs,a)).*kron(wts.',[1 1]);
Tsc = [Tsc; zeros(2*extra,2*npt)];

% crack -> circle displacement
Pk = kron(P,eye(2));
Ucs = (1/c)*(-prob.kd.eval(crk(prob.sm,a),prob.circ) ...
    .*kron(a*a*prob.qm,[1 1]))*Pk;
Ncr = kron(Nself,eye(2));

Rk = eye(size(Ncr,2));
if strcmp(method,'IBPw')
    wr = kron([sqrt(1-xs(:).^2); ones(extra,1)],[1;1]);
    Ncr = wr.*Ncr; Tsc = wr.*Tsc;
end
if strcmp(method,'CP')
    Rk = kron(somegaR(N),eye(2))/(pi^2*Dc);
    Ncr = Ncr*Rk; Ucs = Ucs*Rk;
end

A = [prob.Acc, Ucs; Tsc, Ncr];
rhs = [prob.g/c; zeros(size(Ncr,1),1)];
[x,~,~,iter] = gmres(A,rhs,[],1e-12,size(A,1));
sig = x(1:2*npt); cr = Rk*x(2*npt+1:end);

% errors
wf = reshape(Pk*cr,2,[]);
wexf = prob.wex(prob.sm);
r.errw = max(vecnorm(wf - wexf))/max(vecnorm(wexf));
wt = reshape(kron(Ptip,eye(2))*cr,2,[]);
r.errtip = max(vecnorm(wt - prob.wex([-1 1])))/max(vecnorm(prob.wex([-1 1])));

[xx,yy] = meshgrid(linspace(-0.95,0.95,60));
tg = [xx(:).'; yy(:).']; tg = tg(:,vecnorm(tg) < 0.95);
tg = tg(:,distslit(tg,a) > 0.05);
t = []; t.r = tg;
uex = prob.uex(t);
u = chunkerkerneval(prob.chnkr,prob.kdalt,sig,tg) + ucrack(prob,wf,t);
r.erru = max(vecnorm(reshape(u-uex,2,[])))/max(vecnorm(reshape(uex,2,[])));
r.condA = cond(A);
r.condcr = cond(Ncr);
r.iter = iter(end);
if keep
    r.sig = sig; r.wf = wf; r.prob = prob;
end
end


function plotsol(r,a)
prob = r.prob;
nplot = 120;
xs = linspace(-1,1,nplot);
[xx,yy] = meshgrid(xs,xs);
tp = [xx(:).'; yy(:).'];
in = vecnorm(tp) < 0.99 & abs(tp(2,:)) > 1e-12;
t = []; t.r = tp(:,in);
u = chunkerkerneval(prob.chnkr,prob.kdalt,r.sig,t.r) + ucrack(prob,r.wf,t);
uex = prob.uex(t);
u = reshape(u,2,[]);
err = vecnorm(u - reshape(uex,2,[]))/max(vecnorm(reshape(uex,2,[])));
fplot = {u(1,:), u(2,:), vecnorm(u), log10(err)};
tlabs = {'$u_1$', '$u_2$', '$|u|$', '$\log_{10}$(error)'};
figure(1); clf
for ip = 1:4
    subplot(2,2,ip)
    zz = nan(size(xx));
    zz(in) = fplot{ip};
    h = pcolor(xx,yy,zz);
    set(h,'EdgeColor','none')
    hold on
    plot(prob.chnkr,'k','LineWidth',2)
    plot([-a a],[0 0],'w','LineWidth',2)
    axis equal tight
    set(gca,"box","off","Xtick",[],"Ytick",[]);
    colorbar()
    title(tlabs{ip},'Interpreter','latex','FontSize',12)
end
if exist('sgtitle','file')
    sgtitle(['disk with a traction-free slit, Dirichlet data on the ' ...
        'circle (IBPw, $N = 32$)'],'Interpreter','latex')
end
end


function u = ucrack(prob,wvals,t)
% crack field u_c = -D[phi], phi = a sqrt(1-s^2) w(s), at targets t,
% w given on the fine Gauss-Chebyshev nodes prob.sm (evaluated in blocks)
a = prob.a;
dens = reshape(wvals.*(a*a*prob.qm),[],1);
nt = size(t.r,2); u = zeros(2*nt,1);
src = crk(prob.sm,a);
for i0 = 1:2000:nt
    ii = i0:min(nt,i0+1999);
    tb = []; tb.r = t.r(:,ii);
    u(2*ii(1)-1:2*ii(end)) = -prob.kd.eval(src,tb)*dens;
end
end


function s = crk(x,a)
% points on the slit (-a,0)-(a,0) at parameters x in [-1,1], normal (0,1)
x = x(:).';
s = []; s.r = [a*x; 0*x]; s.n = [0*x; 1+0*x];
end


function d = distslit(t,a)
inside = abs(t(1,:)) <= a;
d = abs(t(2,:)).*inside + ...
    vecnorm(t - [a*sign(t(1,:)); 0*t(1,:)]).*(~inside);
end


function [s,q] = gc2(n)
% Gauss-Chebyshev (2nd kind): int sqrt(1-s^2) f(s) ds ~ sum q_j f(s_j)
th = (1:n)*pi/(n+1); s = cos(th); q = pi/(n+1)*sin(th).^2;
end


function U = umat(s,n)
% U(j,k+1) = U_k(s_j), k = 0..n-1
s = s(:); U = zeros(numel(s),n);
U(:,1) = 1;
if n > 1, U(:,2) = 2*s; end
for k = 3:n, U(:,k) = 2*s.*U(:,k-1) - U(:,k-2); end
end


function S = somegaR(N)
% w at the zeros of U_N from mu at the zeros of T_N, w = S_omega[mu] =
% int log|s-y| mu(y)/sqrt(1-y^2) dy, using
% int log|x-y| T_k(y)/sqrt(1-y^2) dy = -pi/k T_k(x) (k > 0), -pi log 2 (k = 0)
th = (2*(1:N)-1)*pi/(2*N);
Tm = cos((0:N-1).*th(:));
Dm = (2/N)*Tm.'; Dm(1,:) = Dm(1,:)/2;
lamk = [-pi*log(2), -pi./(1:N-1)];
sj = cos((1:N)*pi/(N+1));
Tx = cos((0:N-1).*acos(sj(:)));
S = Tx*diag(lamk)*Dm;
end
