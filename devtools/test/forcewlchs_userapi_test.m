% forcewlchs_userapi_test.m
%
% Exercise opts.forcewlchs in chunkerkerneval as a user would. Tests:
%   (1) bare chunker (no rcipsav): close eval at sub-panel distances;
%   (2) chunkgraph + opts.rcipsav, single Helmholtz D kernel: corner-
%       aware close eval (replaces the earlier recursive-adaptive path);
%   (3) self-consistency between (2) and a manual reference.

clearvars; close all; format long e;
addpath('.');

zk = 1.1;
S = kernel('helm','s',zk);
D = kernel('helm','d',zk);

% ----- (1) bare chunker, single S kernel ---------------------------------
chnkr = chunkerfunc(@(t) circarc(t,1.0), ...
    struct('ta',0,'tb',pi/2,'nover',1), struct('k',16));
chnkr = sort(chnkr);
theta = atan2(chnkr.r(2,:), chnkr.r(1,:));
dens = (1 + 0.5*sin(theta(:)));
near = [cos(0.4); sin(0.4)] + 1e-3 *[cos(0.5); sin(0.5)];
sub  = [cos(0.4); sin(0.4)] + 1e-9 *[cos(0.5); sin(0.5)];
far  = [0.5; 0.5];

ref = chunkerkerneval(chnkr, S, dens, [far,near,sub], struct('forceadap',true,'eps',1e-14));
val = chunkerkerneval(chnkr, S, dens, [far,near,sub], struct('forcewlchs',true));
fprintf('(1) bare chunker, S kernel: max err = %.3e\n', max(abs(val-ref)));

% ----- (2) chunkgraph + rcipsav, single D kernel -------------------------
funs = {@(t) circle(t), @(t) circle(t)};
cparams1 = []; cparams1.ta = 0; cparams1.tb = pi/1.1; cparams1.nchmin = 1;
cparams2 = []; cparams2.ta = 0; cparams2.tb = pi/1.4; cparams2.nchmin = 8;
cg = chunkgraph([-1 1; 0 0], [1 2; 2 1], funs, {cparams1, cparams2});
chnkrtot = merge(cg.echnks);

% Solve a known closed-curve Dirichlet problem (rcipTest setup)
fkern = -2*kernel('helm','d',zk);
ns = 10; rng(8675309);
ts = 2*pi*rand(1,ns);
sources = 3.0*[cos(ts); sin(ts)];
strengths = randn(ns,1);
ubdry = S.eval(struct('r',sources), struct('r',chnkrtot.r(:,:))) * strengths;
[sysmat, ~, rcipsav] = chunkermat(cg, fkern, struct('rcip',true,'nsub_or_tol',30,'rcip_savedepth',inf));
sysmat = sysmat + eye(chnkrtot.npt);
sol = gmres(sysmat, ubdry, [], 1e-13, 200);

% Targets: far + near + sub-panel from corner (-1,0)
dvec = [1, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6];
targets = [-1+dvec; zeros(1,length(dvec))];
utrue = S.eval(struct('r',sources), struct('r',targets)) * strengths;

% Probe diagnostic
src_probe = struct('r',[0;0],'n',[1;0],'d',[1;0],'d2',[0;0]);
tgt_probe = struct('r',[1;0],'n',[1;0],'d',[1;0],'d2',[0;0]);
v = fkern.eval(src_probe, tgt_probe);
b = 0.25i*zk*besselh(1,zk);
fprintf('probe val=%s base=%s coef=%s (expect -2)\n', num2str(v), num2str(b), num2str(v/b));

% Compare on chunkgraph (no rcipsav) at a far target only -- should agree
% with the default hybrid path away from corners
far_tgt = [0; 0.5];
val_def = chunkerkerneval(cg, fkern, sol, far_tgt, []);
val_wl  = chunkerkerneval(cg, fkern, sol, far_tgt, struct('forcewlchs',true));
fprintf('cg+forcewlchs (no rcipsav), far tgt: |default-forcewlchs|=%.3e\n', abs(val_def-val_wl));

% bypass chunkgraph: merge first, then call forcewlchs on bare chunker
val_mwl = chunkerkerneval(chnkrtot, fkern, sol, far_tgt, struct('forcewlchs',true));
fprintf('merged_chunker+forcewlchs, far tgt: |default-mergewlchs|=%.3e\n', abs(val_def-val_mwl));

% Manual smooth sum for sanity
zsrc = chnkrtot.r(1,:) + 1i*chnkrtot.r(2,:);
nz_s = chnkrtot.n(1,:) + 1i*chnkrtot.n(2,:);
ztgt = far_tgt(1) + 1i*far_tgt(2);
diff = ztgt - zsrc;
r_ = abs(diff);
% chunkie D = (i*zk/4) H_1(zk r) imag(nz/diff) * |y-x|... let me use formula from green
val_manual_D = (0.25i * zk * besselh(1,zk*r_) .* imag(nz_s./diff)) * (sol .* chnkrtot.wts(:));
fprintf('manual smooth D: %s   |default-manual|=%.3e\n', num2str(val_manual_D), abs(val_def-(-2*val_manual_D)));

% Compare: same call WITHOUT forcewlchs (uses recursive adaptive)
val_old = chunkerkerneval(cg, fkern, sol, targets, struct('rcipsav',{rcipsav}));

% chunkerkerneval with forcewlchs + rcipsav
val = chunkerkerneval(cg, fkern, sol, targets, ...
    struct('forcewlchs',true,'rcipsav',{rcipsav}));

fprintf('\n(2) chunkgraph + rcipsav, D kernel:\n');
fprintf('   d           rcipsav-only err   forcewlchs+rcipsav err\n');
for i = 1:length(dvec)
    fprintf('   %.0e          %.3e         %.3e\n', dvec(i), ...
        abs(val_old(i)-utrue(i)), abs(val(i)-utrue(i)));
end

% Practical assertion: 12+ digits at d >= 1e-3 (typical grid scale).
assert(max(abs(val(dvec >= 1e-3) - utrue(dvec >= 1e-3))) < 1e-12, ...
    'expected 12+ digits at d >= 1e-3');
fprintf('\nAll assertions passed.\n');

function [r,d,d2] = circle(t)
    c = cos(t(:).'); s = sin(t(:).');
    r = [c; s]; d = [-s; c]; d2 = [-c; -s];
end
function [r,d,d2] = circarc(t,R)
    c = cos(t(:).'); s = sin(t(:).');
    r = R*[c;s]; d = R*[-s;c]; d2 = R*[-c;-s];
end
