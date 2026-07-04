function prm = params(opts)
%CHNK.LNSF2D.PARAMS derived parameters for the 2D time-harmonic linearized
% Navier-Stokes-Fourier (thermoviscous acoustic) system, time convention
% exp(-i omega t):
%
%   mu Lap v + (lam+mu) grad div v + i om rho0 v - eta_v grad T = 0
%   kap Lap T + i om rho0 cv T - eta_T div v = 0
%
% with lam = muB - (2/3) mu + i/(om betaT), eta_v = alpha/betaT,
% eta_T = alpha T0/betaT. The two longitudinal (Kirchhoff) wavenumbers
% ka, kh are the roots of
%
%   s^2 - (K1 + K2 - eta_v eta_T/((lam+2mu) kap)) s + K1 K2 = 0,
%   K1 = i om rho0/(lam+2mu),  K2 = i om rho0 cv/kap,
%
% and kv^2 = i om rho0/mu is the transverse (vortical) wavenumber. All
% wavenumbers are taken with nonnegative imaginary part.
%
% Syntax: prm = chnk.lnsf2d.params(); prm = chnk.lnsf2d.params(opts);
%
% opts fields (with defaults chosen so that the wavenumber ratio is ~5):
%   rho0 (1), om (1), gam (1.4), cv (1), c (1), T0 (1),
%   mu (0.04), muB (0.05), kap (0.0625)
%
% Output struct contains the inputs plus lam, etav, etaT, ka, kh, kv,
% ka2, kh2, kv2, ba, bh, Ba, Bh, cvT, cTv, c0, theta (pseudostress
% parameter), iomrho0, K1, K2, betaT, alpha, cp.
%
% see also CHNK.LNSF2D.KERN

if nargin < 1, opts = []; end
rho0 = getf(opts, 'rho0', 1.0); om = getf(opts, 'om', 1.0);
gam = getf(opts, 'gam', 1.4); cv = getf(opts, 'cv', 1.0);
c = getf(opts, 'c', 1.0); T0 = getf(opts, 'T0', 1.0);
mu = getf(opts, 'mu', 0.04); muB = getf(opts, 'muB', 0.05);
kap = getf(opts, 'kap', 0.0625);

cp = gam*cv;
betaT = gam/(rho0*c*c);
alpha = sqrt((cp - cv)*rho0*betaT/T0);
lam = muB - 2*mu/3 + 1i/(om*betaT);
etav = alpha/betaT;
etaT = alpha*T0/betaT;
K1 = 1i*om*rho0/(lam + 2*mu);
K2 = 1i*om*rho0*cv/kap;
S = K1 + K2 - etav*etaT/((lam+2*mu)*kap);
P = K1*K2;
disc = sqrt(S*S - 4*P);
s1 = (S + disc)/2; s2 = (S - disc)/2;
if abs(s1) < abs(s2), ka2 = s1; kh2 = s2;
else, ka2 = s2; kh2 = s1; end
ka = sqrtu(ka2); kh = sqrtu(kh2); kv = sqrtu(1i*om*rho0/mu);

% ka2 is within O(E) of K1 and kh2 within O(E) of K2, where
% E = etav etaT/((lam+2mu) kap): forming the small differences directly
% loses ~8 digits, which the T rows inherit through Ba and bh. The
% dispersion relation gives exact product forms:
%   (ka2-K1)(kh2-K1) = K1 E,   (kh2-K2)(ka2-K2) = K2 E.
E = etav*etaT/((lam+2*mu)*kap);
dka_K1 = K1*E/(kh2 - K1);         % = ka2 - K1, cancellation-free
dkh_K2 = K2*E/(ka2 - K2);         % = kh2 - K2, cancellation-free
Aa = (ka2 - K2)/(ka2 - kh2); Ah = dkh_K2/(kh2 - ka2);
prm = struct();
prm.rho0 = rho0; prm.om = om; prm.gam = gam; prm.cv = cv; prm.cp = cp;
prm.c = c; prm.T0 = T0; prm.mu = mu; prm.muB = muB; prm.kap = kap;
prm.betaT = betaT; prm.alpha = alpha; prm.lam = lam;
prm.etav = etav; prm.etaT = etaT;
prm.K1 = K1; prm.K2 = K2;
prm.ka = ka; prm.kh = kh; prm.kv = kv;
prm.ka2 = ka2; prm.kh2 = kh2; prm.kv2 = kv*kv;
prm.ba = K1*Aa/ka2; prm.bh = K1*Ah/kh2;
prm.Ba = dka_K1/(ka2 - kh2); prm.Bh = 1.0 - prm.Ba;
prm.cvT = etav/((lam+2*mu)*kap*(ka2 - kh2));
prm.cTv = etaT/((lam+2*mu)*kap*(ka2 - kh2));
prm.iomrho0 = 1i*om*rho0;
prm.c0 = 1/(1i*om*rho0);
prm.theta = mu*(lam + mu)/(lam + 3*mu);

end

function v = getf(o, f, dflt)
if isfield(o, f), v = o.(f); else, v = dflt; end
end

function s = sqrtu(w)
s = sqrt(w);
if imag(s) < 0, s = -s; end
end
