function blk = momentside(prm, S, Itab, Jtab, rs, ds, xi0, sigma, Wside, ...
                          eta, tscale, k, s8, w8)
%CHNK.LNSF2D.MOMENTSIDE moment-integrated singular-piece contribution.
%
% One side (sign sigma) of the singular piece around the self target at
% node parameter xi0: the integral over xi in [xi0, xi0 + sigma Wside]
% of (D + eta S)(x(xi0), y(xi)) * (density basis) * |c'|. The kernel's
% radial profiles (log-series from CHNK.LNSF2D.PROFSERIES) are
% contracted against exact moments (CHNK.LNSF2D.MOMTABS) of a Legendre
% fit of the geometry in the scaled radial variable v = r/R. The
% kernel's large collar values are never sampled; the only conditioning
% cost is e^{|k| R} of the series (keep |k_max| R <= ~7).
%
% Output blk is (3 x 3k), density interleaved (v1,v2,T) per node, with
% the tscale conjugation applied (T-row / tscale, T-columns * tscale).
%
% see also CHNK.LNSF2D.NEARCORRECT, CHNK.LNSF2D.PROFSERIES

PF = 22;                                   % fit degree: high enough to
                                           % resolve the degree-15 density
                                           % basis through the u-warp even
                                           % when the piece spans the chunk
[xg, ~] = lege.exps(PF+1);
vs = 0.5*(xg(:).' + 1);                    % fit samples in (0,1)
persistent VI PFsave
if isempty(VI) || PFsave ~= PF
    pol = lege.pols(xg(:), PF);
    if size(pol, 1) ~= PF+1, pol = pol.'; end
    VI = inv(pol.');                       % values at samples -> coefs
    PFsave = PF;
end

% ---- radius at the outer end, then Newton for w(v R) ----
uofw = @(w) uz(w, sigma, xi0, rs, ds, s8, w8, k);
[R, ~] = uofw(Wside); R = R(1);
w = vs*Wside;
for it = 1:4
    [u, z] = uofw(w);
    P1 = lege.matrin(k, xi0 + sigma*w(:));
    dfin = (P1*ds.').';
    up = -sigma*(z(1,:).*dfin(1,:) + z(2,:).*dfin(2,:))./u;   % du/dw
    w = w - (u - vs*R)./up;
    w = min(max(w, 1e-300), Wside);
end
[u, z] = uofw(w);
taus = xi0 + sigma*w;
Pm = lege.matrin(k, taus(:));
dfin = (Pm*ds.').';
sp = sqrt(sum(dfin.^2, 1));
nn = [dfin(2,:); -dfin(1,:)]./sp;
zh = z./u;
et = stablee(w, sigma, xi0, ds, s8, w8, k)./(sp.*u.^2);  % e/u
up = -sigma*(z(1,:).*dfin(1,:) + z(2,:).*dfin(2,:))./u;
jac = sp./up;                              % sp * |dtau/du| (both > 0)
PHI = Pm;                                  % (PF+1, k) density basis

mu = prm.mu; th = prm.theta; lmt = prm.lam + mu - th;
blk = zeros(3, 3*k);

one = ones(size(u));
for i = 1:2
    for j = 1:2
        if i == j
            blk = accum(blk, i, j, S.ua, -et, jac, PHI, VI, Itab, Jtab, R, PF);
            blk = accum(blk, i, j, S.gP, eta*one, jac, PHI, VI, Itab, Jtab, R, PF);
        end
        blk = accum(blk, i, j, S.b,  -nn(i,:).*zh(j,:), jac, PHI, VI, Itab, Jtab, R, PF);
        blk = accum(blk, i, j, S.c,  -nn(j,:).*zh(i,:), jac, PHI, VI, Itab, Jtab, R, PF);
        blk = accum(blk, i, j, S.ud, -zh(i,:).*zh(j,:).*et, jac, PHI, VI, Itab, Jtab, R, PF);
        blk = accum(blk, i, j, S.Q1, eta*zh(i,:).*zh(j,:), jac, PHI, VI, Itab, Jtab, R, PF);
    end
    blk = accum(blk, i, 3, S.Pd,   prm.kap*prm.cvT*nn(i,:), jac, PHI, VI, Itab, Jtab, R, PF);
    blk = accum(blk, i, 3, S.uQd,  prm.kap*prm.cvT*zh(i,:).*et, jac, PHI, VI, Itab, Jtab, R, PF);
    blk = accum(blk, i, 3, S.chip, -eta*prm.cvT*zh(i,:), jac, PHI, VI, Itab, Jtab, R, PF);
end
for j = 1:2
    blk = accum(blk, 3, j, S.Pd,   prm.cTv*((mu+th) + 2*lmt)*nn(j,:), jac, PHI, VI, Itab, Jtab, R, PF);
    blk = accum(blk, 3, j, S.uQd,  prm.cTv*(mu+th)*zh(j,:).*et, jac, PHI, VI, Itab, Jtab, R, PF);
    blk = accum(blk, 3, j, S.Qd,   prm.cTv*lmt*nn(j,:), jac, PHI, VI, Itab, Jtab, R, PF);
    blk = accum(blk, 3, j, S.chip, -eta*prm.cTv*zh(j,:), jac, PHI, VI, Itab, Jtab, R, PF);
end
blk = accum(blk, 3, 3, S.upsip, -et, jac, PHI, VI, Itab, Jtab, R, PF);
blk = accum(blk, 3, 3, S.psi,   (eta/prm.kap)*one, jac, PHI, VI, Itab, Jtab, R, PF);

if tscale ~= 1.0
    blk(3, :) = blk(3, :)/tscale;
    blk(:, 3:3:end) = blk(:, 3:3:end)*tscale;
end
end

function blk = accum(blk, row, col0, prof, struct_, jac, PHI, VI, Itab, Jtab, R, PF)
g = (struct_.*jac).'.*PHI;                 % (PF+1, k)
gh = VI*g;                                 % Legendre coefficients of g(v)
q = prof.o + (0:numel(prof.P)-1).';
Rp = R.^q;
Pt = (prof.P + prof.Q*log(R)).*Rp;
Qt = prof.Q.*Rp;
qq = min(q, size(Itab,1)-1);
row16 = (Pt.'*Itab(qq+1, 1:PF+1) + Qt.'*Jtab(qq+1, 1:PF+1));
blk(row, col0:3:end) = blk(row, col0:3:end) + R*(row16*gh);
end

function [u, z] = uz(w, sigma, xi0, rs, ds, s8, w8, k) %#ok<INUSL>
% stable displacement z = x - y and radius u at |dxi| = w along side sigma
w = w(:).';
dxi = sigma*w;
taus = xi0 + dxi(:)*s8;                    % (nw, 8)
P8 = lege.matrin(k, taus(:));
dtau = (P8*ds.').';
dtau = reshape(dtau, 2, numel(w), numel(s8));
zint = zeros(2, numel(w));
for m = 1:numel(w8)
    zint = zint + w8(m)*dtau(:,:,m);
end
z = -dxi.*zint;
u = sqrt(sum(z.^2, 1));
end

function ndz = stablee(w, sigma, xi0, ds, s8, w8, k)
% n(y) . z in cancellation-free form:
%   n.z = dxi^2/sp * int_0^1 int_0^1 (s-1) c'(tau) x c''(tau + t(s-1)dxi)
% (c' x c' = 0 removes the eps-size absolute error of the plain dot at
% small r; c'' is exact for polynomial chunks). Returns sp * e * u, i.e.
% divide by (sp u^2) for e/u.
w = w(:).';
dxi = sigma*w;
% Legendre coefficients of c' and their derivative (c'')
[~, ~, uu, ~] = lege.exps(k);
cds = (uu*ds.').';                         % (2, k) coefficients of c'
cdd = zeros(2, k);
for n = k-2:-1:0                           % derivative of Leg. expansion
    hi = n+1:2:k-1;                        % j > n with j - n odd
    cdd(:, n+1) = (2*n+1)*sum(cds(:, hi+1), 2);
end
nw = numel(w); n8 = numel(s8);
ss = s8(:).';                               % quadrature on [0,1]
poly = lege.pols(xi0 + dxi(:), k-1);
if size(poly, 1) ~= k, poly = poly.'; end
cpv = cds*poly;                             % (2, nw) c' at y
acc = zeros(1, nw);
for m1 = 1:n8                              % outer: s
    sarg = ss(m1);
    inner = zeros(1, nw);
    for m2 = 1:n8
        targ = (xi0 + dxi(:)) + s8(m2)*(sarg - 1)*dxi(:);   % (nw,1)
        pol = lege.pols(targ, k-1);
        if size(pol, 1) ~= k, pol = pol.'; end
        cppv = cdd*pol;                    % (2, nw) c'' at targ
        cr = cpv(1,:).*cppv(2,:) - cpv(2,:).*cppv(1,:);
        inner = inner + w8(m2)*cr;
    end
    acc = acc + w8(m1)*(sarg - 1)*inner;
end
ndz = (dxi.^2).*acc;                       % = sp * (n . z)
end
