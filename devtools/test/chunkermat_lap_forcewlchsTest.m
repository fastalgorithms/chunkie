function chunkermat_lap_forcewlchsTest()
%CHUNKERMAT_LAP_FORCEWLCHSTEST  Verify chunkermat opts.forcewlchs dispatches
%   correctly to chnk.kernsplit Laplace self/adj corrections (s, d, sp, dp).
%
%   The test checks operator action against analytic identities on a unit
%   circle.  We don't compare matrix entries against chunkermat default
%   because the two quadrature schemes (chunkermat default's smooth-GL +
%   splitinfo vs forcewlchs's smooth-GL + Helsing log moments) encode the
%   same operator differently and need not agree entry-wise.

cparams = []; cparams.ta = 0; cparams.tb = 2*pi; cparams.ifclosed = true;
cparams.maxchunklen = 2*pi/8;
pref = []; pref.k = 16;
chnkr = chunkerfunc(@(t) circle_param(t), cparams, pref);
N = chnkr.npt;
sig1 = ones(N, 1);

% --- 's' single-layer ---
% On unit circle, S[1] = -1/(2*pi) * integral log(|x-y|) dy = 0.
S_wl = chunkermat(chnkr, kernel('lap','s'), ...
                  struct('quad','native','forcewlchs',true));
val_S = S_wl * sig1;
err_S = norm(val_S, 'inf');
fprintf('  lap s  : ||S*1||_inf  = %.3e (expect 0 on unit circle)\n', err_S);
assert(err_S < 1e-12);

% --- 'd' double-layer: D[1] = -1/2 (interior limit on closed smooth curve) ---
D_wl = chunkermat(chnkr, kernel('lap','d'), ...
                  struct('quad','native','forcewlchs',true));
val_D = D_wl * sig1;
err_D = norm(val_D - (-0.5), 'inf');
fprintf('  lap d  : ||D*1 - (-1/2)||_inf = %.3e\n', err_D);
assert(err_D < 1e-12);

% --- 'sp' normal derivative of single layer: K*[1] = -1/2 in chunkie convention.
%     SLP(1)|interior = 0 (harmonic, zero on boundary).  Interior normal-deriv
%     trace = +1/2 + K*[1] = 0 (Colton-Kress jump rel.) => K*[1] = -1/2.
SP_wl = chunkermat(chnkr, kernel('lap','sp'), ...
                   struct('quad','native','forcewlchs',true));
val_SP = SP_wl * sig1;
err_SP = norm(val_SP - (-0.5), 'inf');
fprintf('  lap sp : ||S''*1 - (-1/2)||_inf = %.3e\n', err_SP);
assert(err_SP < 1e-11);

% --- 'dp' hypersingular: D'[1] = 0 (since D maps constants to constants) ---
DP_wl = chunkermat(chnkr, kernel('lap','dp'), ...
                   struct('quad','native','forcewlchs',true));
val_DP = DP_wl * sig1;
err_DP = norm(val_DP, 'inf');
fprintf('  lap dp : ||D''*1||_inf = %.3e\n', err_DP);
assert(err_DP < 1e-10);

% --- Boolean form on a scalar-multiplied kernel: 3.7 * S[1] = 0 ---
kern_scaled = 3.7 * kernel('lap','s');
S_scaled = chunkermat(chnkr, kern_scaled, ...
                      struct('quad','native','forcewlchs',true));
err_scaled = norm(S_scaled * sig1, 'inf');
fprintf('  3.7*s  : ||3.7*S*1||_inf = %.3e\n', err_scaled);
assert(err_scaled < 1e-11);

fprintf('PASS: chunkermat forcewlchs lap dispatch verified for s/d/sp/dp\n');
end

function [r, d, d2] = circle_param(t)
ct = cos(t); st = sin(t);
r  = [ct(:).';  st(:).'];
d  = [-st(:).'; ct(:).'];
d2 = [-ct(:).'; -st(:).'];
end
