function chunkermat_stok_forcewlchsTest()
%CHUNKERMAT_STOK_FORCEWLCHSTEST  Verify chunkermat opts.forcewlchs
%   dispatches to chnk.kernsplit Stokes self/adj corrections for both
%   the Stokes velocity SLP and DLP kernels.
%
%   Test: builds the forcewlchs Stokes matrices and compares operator
%   action against chunkermat's default scheme on a smooth density
%   (avoiding the radial-normal direction which is effectively in the
%   Stokes SLP null space on a unit circle).

cparams = []; cparams.ta = 0; cparams.tb = 2*pi; cparams.ifclosed = true;
cparams.maxchunklen = 2*pi/8;
pref = []; pref.k = 16;
chnkr = chunkerfunc(@(t) circle_param(t), cparams, pref);
mu = 1.7;

N = chnkr.npt;
r  = chnkr.r(:,:);
ang = atan2(r(2,:), r(1,:)).';
sig = zeros(2*N, 1);
sig(1:2:end) = cos(2*ang);
sig(2:2:end) = sin(2*ang);

% --- SLP (svel) ---
S_def = chunkermat(chnkr, kernel('stok','svel',mu));
S_wl  = chunkermat(chnkr, kernel('stok','svel',mu), ...
                   struct('quad','native','forcewlchs',true));
err = norm(S_def*sig - S_wl*sig) / norm(S_def*sig);
fro_err = norm(S_def - S_wl, 'fro') / norm(S_def, 'fro');
fprintf('  svel: action %.3e (rel),  Frob %.3e (rel)\n', err, fro_err);
assert(err < 1e-9);
assert(fro_err < 1e-9);

% --- DLP (dvel): bounded kernel; default already smooth-GL.  forcewlchs
%     reproduces it via analytic diagonal limit. ---
D_def = chunkermat(chnkr, kernel('stok','dvel',mu));
D_wl  = chunkermat(chnkr, kernel('stok','dvel',mu), ...
                   struct('quad','native','forcewlchs',true));
% Use a non-null-space density: constant (1,0).  The Stokes DLP interior
% trace identity gives D*sigma_const = -sigma_const/2.
sig_const = zeros(2*N,1); sig_const(1:2:end) = 1;
errD = norm(D_def*sig_const - D_wl*sig_const) / norm(D_def*sig_const);
fro_errD = norm(D_def - D_wl, 'fro') / norm(D_def, 'fro');
fprintf('  dvel: action %.3e (rel),  Frob %.3e (rel)\n', errD, fro_errD);
assert(errD < 1e-9);
assert(fro_errD < 1e-9);

% --- D*1 = -1/2 identity (interior trace of Stokes DLP). ---
sig1 = ones(2*N, 1);
val_D1 = D_wl * sig1;
err_id = norm(val_D1 + 0.5*sig1, 'inf');
fprintf('  dvel: ||D*1 + 1/2||_inf = %.3e\n', err_id);
assert(err_id < 1e-12);

% --- Boolean form on a scalar-multiplied Stokes svel ---
kern_scaled = 2.5 * kernel('stok','svel',mu);
S_scaled_def = chunkermat(chnkr, kern_scaled);
S_scaled_wl  = chunkermat(chnkr, kern_scaled, ...
                          struct('quad','native','forcewlchs',true));
err2 = norm(S_scaled_def * sig - S_scaled_wl * sig) / norm(S_scaled_def * sig);
fprintf('  2.5*svel scaled match: %.3e (rel)\n', err2);
assert(err2 < 1e-9);

fprintf('PASS: chunkermat forcewlchs stok dispatch works for svel and dvel\n');
end

function [r, d, d2] = circle_param(t)
ct = cos(t); st = sin(t);
r  = [ct(:).';  st(:).'];
d  = [-st(:).'; ct(:).'];
d2 = [-ct(:).'; -st(:).'];
end
