chunkermat_lnsf2dTest0();

function chunkermat_lnsf2dTest0()
%CHUNKERMAT_LNSF2DTEST
%
% Tests for the 2D time-harmonic linearized Navier-Stokes-Fourier
% (thermoviscous acoustic) kernels in chnk.lnsf2d:
%
%   (1) reference-value check of the S and pseudostress-D kernels against
%       an independently verified implementation
%   (2) finite-difference check that kernel columns satisfy the PDE
%   (3) method-of-manufactured-solutions test of the exterior Dirichlet
%       (no-slip + isothermal) combined-field BIE
%          (1/2 I + D + eta S) psi = -u^in on Gamma,   u = (D + eta S) psi
%       with manufactured solution given by interior point sources
%   (4) accuracy as a function of the chunk length measured in units of
%       the boundary-layer decay length delta = 1/Im(kh)
%
% The default parameters give |kv|/|ka| ~ 5, |kh|/|ka| ~ 4.8.

iseed = 8675309;
rng(iseed);

prm = chnk.lnsf2d.params();
fprintf('ka = %8.5f%+8.5fi\n', real(prm.ka), imag(prm.ka));
fprintf('kh = %8.5f%+8.5fi\n', real(prm.kh), imag(prm.kh));
fprintf('kv = %8.5f%+8.5fi\n', real(prm.kv), imag(prm.kv));

% dispersion relation residuals
resa = prm.ka2^2 - (prm.K1+prm.K2-prm.etav*prm.etaT/((prm.lam+2*prm.mu)*prm.kap))*prm.ka2 + prm.K1*prm.K2;
resh = prm.kh2^2 - (prm.K1+prm.K2-prm.etav*prm.etaT/((prm.lam+2*prm.mu)*prm.kap))*prm.kh2 + prm.K1*prm.K2;
assert(abs(resa) < 1e-12 && abs(resh) < 1e-10);
assert(abs(prm.ba + prm.bh - 1) < 1e-13 && abs(prm.Ba + prm.Bh - 1) < 1e-13);

% ---------------------------------------------------------------
% (1) reference values (independently verified implementation; the same
% kernels were validated there by PDE finite differences, jump-relation
% tests, and an MMS solve converging to ~1e-7)
% ---------------------------------------------------------------
srcinfo = []; srcinfo.r = [0.05; -0.1]; srcinfo.n = [cos(0.3); sin(0.3)];
t1 = []; t1.r = [0.6; 0.45];
t2 = []; t2.r = [0.05 + 0.01*cos(1.1); -0.1 + 0.01*sin(1.1)];

S_t1 = [0.0223962010367964+0.000872492338982123i, 0.0882479905973797+0.339065982411661i, 0.0975831597915792+0.0302254966489512i; 0.0882479905973797+0.339065982411661i, 0.0223962010367964+0.000872492338982234i, 0.0975831597915793+0.0302254966489512i; 0.0975831597915792+0.0302254966489512i, 0.0975831597915793+0.0302254966489512i, -0.0567128611214756+0.0203802998790955i];
D_t1 = [0.084762774674118+0.0344377063749013i, 0.0547764087900867+0.0675389427242625i, 0.00308246775155403-0.00658913278139032i; 0.139748151420657+0.0994728091536361i, 0.0389369665322228+0.0494235153055743i, 0.0103991674930082-0.00432285155019179i; 0.13618773182063-0.0160034322861888i, 0.0542450479665435-0.00719227684946963i, -0.0194507015728641+0.0223498147849966i];
S_t2 = [5.81970364452604+1.07689660405842i, 0.797421492270555+0.0460479921974514i, 0.00577998079379119-0.0148582276362155i; 0.797421492270555+0.0460479921974514i, 6.9805831163644+1.14393288208647i, 0.011356273083313-0.0291928462378532i; 0.00577998079379119-0.0148582276362155i, 0.011356273083313-0.0291928462378532i, 8.10993021974265+1.89931935648685i];
D_t2 = [4.69741667394372-0.671607387594815i, 8.76116868197924+0.992080321412984i, -0.0733775857550589+0.177163895839969i; 8.77488585165342+0.956522536109027i, 17.4673282809581+0.735378431478273i, -0.0182180539776239+0.024309595742829i; -3.97184638636358-1.28648861636017i, -1.22086206453966-0.436259044565692i, 11.0770603291503+0.0448829889845465i];

errs = zeros(4,1);
errs(1) = norm(chnk.lnsf2d.kern(prm, srcinfo, t1, 's') - S_t1)/norm(S_t1);
errs(2) = norm(chnk.lnsf2d.kern(prm, srcinfo, t1, 'd') - D_t1)/norm(D_t1);
errs(3) = norm(chnk.lnsf2d.kern(prm, srcinfo, t2, 's') - S_t2)/norm(S_t2);
errs(4) = norm(chnk.lnsf2d.kern(prm, srcinfo, t2, 'd') - D_t2)/norm(D_t2);
fprintf('reference-value rel errs: %.2e %.2e %.2e %.2e\n', errs);
assert(all(errs < 1e-10));

% ---------------------------------------------------------------
% (2) finite-difference PDE check on kernel columns
% ---------------------------------------------------------------
for tp = 'sd'
    for jcol = 1:3
        ucol = @(x) col_of_kern(prm, srcinfo, x, tp, jcol);
        for xt = [[0.6; 0.45], [-0.55; -0.35]]
            [res, u0] = pde_residual(prm, ucol, xt, 4e-3);
            rel = max(abs(res))/max(max(abs(u0)), 1e-12);
            assert(rel < 2e-5);
        end
    end
end
fprintf('FD PDE tests on S and D columns passed\n');

% ---------------------------------------------------------------
% (3),(4) MMS: exterior Dirichlet solve, chunk-length study
% ---------------------------------------------------------------
delta = 1/imag(prm.kh);       % thermal decay length ~ viscous one here
fprintf('decay length delta = %5.3f\n', delta);

% manufactured solution: interior point sources
srcpts = []; srcpts.r = [0.10, -0.20; -0.35, 0.15];
strengths = randn(6,1) + 1i*randn(6,1);

% exterior targets
tt = 2*pi*(0:7)/8;
targs = [2.2*cos(tt), 1.62*cos(tt(1:2:end)+0.3); ...
         2.2*sin(tt), 1.62*sin(tt(1:2:end)+0.3)];
uex = @(pts) chnk.lnsf2d.kern(prm, srcpts, struct('r', pts), 's')*strengths;
utrue = uex(targs);

eta = 1.0;
fkern = kernel(@(s,t) chnk.lnsf2d.kern(prm, s, t, 'c', eta));
fkern.opdims = [3 3];
fkern.sing = 'log';

% fixed uniform chunk counts so that coarse chunks are actually coarse;
% nch = 24 gives max chunk length ~ delta, nch = 4 gives ~ 6 delta
narms = 3; amp = 0.25;
nchs = [24, 16, 12, 8, 6, 4];
relerrs = zeros(size(nchs));
for imc = 1:length(nchs)
    chnkr = chunkerfuncuni(@(t) starfish(t, narms, amp), nchs(imc));
    maxlen = max(chunklen(chnkr));
    sysmat = chunkermat(chnkr, fkern);
    sys = 0.5*eye(size(sysmat,1)) + sysmat;
    bdrypts = reshape(chnkr.r, 2, chnkr.k*chnkr.nch);
    rhs = uex(bdrypts);
    psi = sys \ rhs;
    opts = []; opts.usesmooth = false;
    usol = chunkerkerneval(chnkr, fkern, psi, targs, opts);
    relerrs(imc) = max(abs(usol - utrue))/max(abs(utrue));
    fprintf('nch = %3d  maxchunklen/delta = %5.2f  relerr = %.3e\n', ...
        nchs(imc), maxlen/delta, relerrs(imc));
end
assert(relerrs(1) < 1e-6);
assert(all(relerrs(1:4) < 1e-5));

end

% ------------------------------------------------------------------
function u = col_of_kern(prm, srcinfo, x, tp, jcol)
targ = []; targ.r = x(:);
K = chnk.lnsf2d.kern(prm, srcinfo, targ, tp);
u = K(:, jcol);
end

function [res, u0] = pde_residual(prm, ufun, x, h)
% 4th-order finite-difference application of the PDE operator
c1 = [1, -8, 0, 8, -1]/12;
c2 = [-1, 16, -30, 16, -1]/12;
offs = -2:2;
d1 = @(dim, comp) stencil1(ufun, x, h, dim, comp, c1, offs)/h;
d2f = @(dim, comp) stencil1(ufun, x, h, dim, comp, c2, offs)/h^2;
d11 = @(comp) stencil2(ufun, x, h, comp, c1, offs)/h^2;
u0 = ufun(x);
lap = zeros(3,1); gT = zeros(2,1);
for cc = 1:3, lap(cc) = d2f(1,cc) + d2f(2,cc); end
divv = d1(1,1) + d1(2,2);
graddiv = [d2f(1,1) + d11(2); d11(1) + d2f(2,2)];
gT(1) = d1(1,3); gT(2) = d1(2,3);
res = zeros(3,1);
res(1) = prm.mu*lap(1) + (prm.lam+prm.mu)*graddiv(1) + prm.iomrho0*u0(1) - prm.etav*gT(1);
res(2) = prm.mu*lap(2) + (prm.lam+prm.mu)*graddiv(2) + prm.iomrho0*u0(2) - prm.etav*gT(2);
res(3) = prm.kap*lap(3) + 1i*prm.om*prm.rho0*prm.cv*u0(3) - prm.etaT*divv;
end

function s = stencil1(ufun, x, h, dim, comp, cf, offs)
s = 0;
for m = 1:length(cf)
    if cf(m) == 0, continue; end
    xx = x; xx(dim) = xx(dim) + offs(m)*h;
    u = ufun(xx); s = s + cf(m)*u(comp);
end
end

function s = stencil2(ufun, x, h, comp, cf, offs)
s = 0;
for m1 = 1:length(cf)
    if cf(m1) == 0, continue; end
    for m2 = 1:length(cf)
        if cf(m2) == 0, continue; end
        xx = x; xx(1) = xx(1) + offs(m1)*h; xx(2) = xx(2) + offs(m2)*h;
        u = ufun(xx); s = s + cf(m1)*cf(m2)*u(comp);
    end
end
end
