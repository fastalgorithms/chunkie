function S = profseries(prm)
%CHNK.LNSF2D.PROFSERIES radial log-series of the combined-kernel profiles.
%
% Each profile is sum_j (P(j) + Q(j) log r) r^(o+j-1), stored as a struct
% with fields o (offset), P, Q (column vectors). The mode combinations
% and the pseudostress cancellations are done at coefficient level, so
% no term is larger than e^{|k| r} times the result. Profiles (z = x-y,
% r = |z|, zh = z/r, e = n(y).zh, et = e/r):
%
%   Wvv_ij = -[ ua*et del_ij + b n_i zh_j + c n_j zh_i + ud zh_i zh_j et ]
%   WvT_i  = kap cvT (Pd n_i + uQd zh_i et)
%   WTv_j  = cTv [ (mu+th)(Pd n_j + uQd zh_j et) + lmt (2Pd + Qd) n_j ]
%   WTT    = -upsip et
%   Gvv_ij = gP del_ij + Q1 zh_i zh_j;  GvT_i = -cvT chip zh_i
%   GTv_j  = -cTv chip zh_j;            GTT = psi/kap
%
% (ua = r*a etc.: profiles that multiply e carry a factor r folded in so
% that all series start at nonnegative powers.)
%
% see also CHNK.LNSF2D.MOMENTSIDE, CHNK.LNSF2D.MOMTABS

MT = 48;                       % highest power kept
mu = prm.mu; th = prm.theta; lam = prm.lam;
lmt = lam + mu - th;

Ga = gser(prm.ka, MT); Gh = gser(prm.kh, MT); Gv = gser(prm.kv, MT);

Phi = sadd(sscale(Gv, prm.c0), sscale(Ga, -prm.c0*prm.ba), ...
           sscale(Gh, -prm.c0*prm.bh));
Phi.Q(1) = 0;                  % exact log cancellation (ba + bh = 1)
gam = sscale(Gv, 1/mu);
chi = sadd(Ga, sscale(Gh, -1));
chi.Q(1) = 0;                  % exact log cancellation
psi = sadd(sscale(Ga, prm.Ba), sscale(Gh, prm.Bh));

Phip = sder(Phi); Phipp = sder(Phip); Phippp = sder(Phipp);
gamp = sder(gam);
chip = sder(chi); chipp = sder(chip);
psip = sder(psi);

P1 = sdivr(Phip);
Q1 = sadd(Phipp, sscale(P1, -1));
R  = sdivr(Q1);
Sw = sadd(Phippp, sscale(sdivr(Phipp), -3), sscale(sdivr(sdivr(Phip)), 3));
% gamp + (Lap Phi)' = c0(ba ka^2 Ga' + bh kh^2 Gh') exactly (PDE
% identity); the assembled form cancels two 1/mu-scale series
gpl = sadd(sscale(sder(Ga), prm.c0*prm.ba*prm.ka^2), ...
           sscale(sder(Gh), prm.c0*prm.bh*prm.kh^2));
Pd = sdivr(chip);
Qd = sadd(chipp, sscale(Pd, -1));

a = sadd(sscale(sadd(gamp, R), mu), sscale(R, th));
b = sadd(sscale(R, mu), sscale(sadd(gamp, R), th));
c = sadd(sscale(R, mu+th), sscale(gpl, lmt));
d = sscale(Sw, mu+th);
% the 1/r parts of b and c cancel exactly by the pseudostress identity;
% in floating point the 1/(4 pi mu)-scale coefficients leave an eps/mu
% residue: zero the r^{-1} coefficients exactly
b = zeropow(b, -1); c = zeropow(c, -1);

S = [];
S.ua   = strim(sshift(a));
S.b    = strim(b);
S.c    = strim(c);
S.ud   = strim(sshift(d));
S.gP   = strim(sadd(gam, P1));
S.Q1   = strim(Q1);
S.Pd   = strim(Pd);
S.Qd   = strim(Qd);
S.uQd  = strim(sshift(Qd));
S.chip = strim(chip);
S.upsip= strim(sshift(psip));
S.psi  = strim(psi);
end

function f = gser(k, MT)
% (i/4) H0^(1)(k r) as a log-series (even powers)
n = floor(MT/2) + 1;
cm = zeros(n, 1); cm(1) = 1;
for i = 1:n-1
    cm(i+1) = cm(i)*(-(k/2)^2)/i^2;
end
Hm = [0; cumsum(1./(1:n-1).')];
gamma_e = 0.57721566490153286;
av = (1i/4)*cm - (log(k/2) + gamma_e)*cm/(2*pi) + Hm.*cm/(2*pi);
bv = -cm/(2*pi);
f = []; f.o = 0;
f.P = zeros(MT+1, 1); f.Q = zeros(MT+1, 1);
f.P(1:2:end) = av; f.Q(1:2:end) = bv;
end

function f = sadd(varargin)
o = min(cellfun(@(g) g.o, varargin));
n = max(cellfun(@(g) g.o + numel(g.P), varargin)) - o;
f = []; f.o = o; f.P = zeros(n, 1); f.Q = zeros(n, 1);
for m = 1:nargin
    g = varargin{m};
    i0 = g.o - o;
    f.P(i0+1:i0+numel(g.P)) = f.P(i0+1:i0+numel(g.P)) + g.P;
    f.Q(i0+1:i0+numel(g.Q)) = f.Q(i0+1:i0+numel(g.Q)) + g.Q;
end
end

function f = sscale(g, s)
f = g; f.P = s*g.P; f.Q = s*g.Q;
end

function f = sder(g)
j = (g.o + (0:numel(g.P)-1)).';
f = []; f.o = g.o - 1; f.P = j.*g.P + g.Q; f.Q = j.*g.Q;
end

function f = sdivr(g)
f = g; f.o = g.o - 1;
end

function f = sshift(g)
f = g; f.o = g.o + 1;
end

function f = zeropow(g, pw)
f = g;
idx = pw - g.o + 1;
if idx >= 1 && idx <= numel(g.P)
    f.P(idx) = 0; f.Q(idx) = 0;
end
end

function f = strim(g)
nz = find(g.P ~= 0 | g.Q ~= 0, 1, 'first');
if isempty(nz), nz = 1; end
f = []; f.o = g.o + nz - 1; f.P = g.P(nz:end); f.Q = g.Q(nz:end);
if f.o < 0
    error('CHNK.LNSF2D.PROFSERIES: negative power survived trimming');
end
end
