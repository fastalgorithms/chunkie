% chunkerkerneval_rcipsav_Test.m
%
% Tests the opts.rcipsav code path in chunkerkerneval. Uses the same setup
% as devtools/test/rcipTest.m -- closed curve from two arcs of unit
% circle meeting at two corners, single Helmholtz D kernel BIE -- so the
% analytic field is computable from point sources outside.
%
% Compares:
%   - naive chunkerkerneval (RCIP-unaware): expected to lose accuracy
%     near each corner because dens(starind) is the RCIP-transformed
%     density, only valid for far-field smooth quadrature.
%   - chunkerkerneval with opts.rcipsav: expected to maintain ~1e-12
%     accuracy down to ~1e-5 from a corner (validated in step-2 smoke
%     test using identical inner machinery).

clearvars; close all; format long e;

zk = 1.1;
funs = {@(t) circle(t), @(t) circle(t)};
cparams1 = []; cparams1.ta = 0; cparams1.tb = pi/1.1; cparams1.nchmin = 1;
cparams2 = []; cparams2.ta = 0; cparams2.tb = pi/1.4; cparams2.nchmin = 8;
cg = chunkgraph([-1 1; 0 0], [1 2; 2 1], funs, {cparams1, cparams2});
chnkrtotal = merge(cg.echnks);

fkern = -2*kernel('helm','d',zk);
kerns = kernel('helm','s',zk);

ns = 10; rng(8675309);
ts = 2*pi*rand(1,ns);
sources = 3.0*[cos(ts); sin(ts)];
strengths = randn(ns,1);

ubdry = kerns.eval(struct('r',sources), struct('r',chnkrtotal.r(:,:))) * strengths;

opts = []; opts.rcip = true; opts.nsub_or_tol = 30;
opts.rcip_savedepth = inf;
fprintf('Building system matrix...\n');
[sysmat, ~, rcipsav] = chunkermat(cg, fkern, opts);
sysmat = sysmat + eye(chnkrtotal.npt);

fprintf('GMRES solve...\n');
sol = gmres(sysmat, ubdry, [], 1e-13, 200);

% Targets: a few far ones plus a sweep approaching corner at (-1,0)
ts = 2*pi*rand(1,3);
targets_far = 0.2*[cos(ts); sin(ts)];
targets_far(:,1) = [0; 0.36];
targets_far(:,2) = [-0.95; 0];
targets_far(:,3) = [0.5; 0.2];

dvec = 10.^(-(1:0.5:6));
targets_near = [-1; 0] + [cos(pi/4)*dvec; sin(pi/4)*dvec];

targets = [targets_far, targets_near];
nt = size(targets,2);

% Reference: free-space single-layer evaluation
utarg = kerns.eval(struct('r',sources), struct('r',targets)) * strengths;

% (1) naive: chunkerkerneval on merged chunker, no rcipsav
opts_eval = [];
unum_naive = chunkerkerneval(chnkrtotal, fkern, sol, targets, opts_eval);

% (2) rcipsav-aware: chunkerkerneval on chunkgraph, with rcipsav
opts_eval_rcip = []; opts_eval_rcip.rcipsav = rcipsav;
unum_rcip = chunkerkerneval(cg, fkern, sol, targets, opts_eval_rcip);

err_naive = abs(unum_naive - utarg);
err_rcip  = abs(unum_rcip  - utarg);

fprintf('\nTargets at fixed positions (far from corners):\n');
fprintf('  %-22s %-14s %-14s\n', 'target', 'naive err', 'rcip-aware err');
for i = 1:size(targets_far,2)
    fprintf('  (%+.4f,%+.4f)      %-14.3e %-14.3e\n', ...
        targets_far(1,i), targets_far(2,i), err_naive(i), err_rcip(i));
end
fprintf('\nTargets approaching corner at (-1,0):\n');
fprintf('  %-12s %-14s %-14s\n', 'd from vert', 'naive err', 'rcip-aware err');
for j = 1:size(targets_near,2)
    i = size(targets_far,2) + j;
    fprintf('  %-12.2e %-14.3e %-14.3e\n', dvec(j), err_naive(i), err_rcip(i));
end

% Far targets should be 1e-12 or better.
assert(max(err_rcip(1:size(targets_far,2))) < 1e-12, ...
    'rcip-aware should match naive at far targets');
% Near-corner: 1e-10 down to d=3e-5; degrades gradually below that as
% targets enter the sub-deepest-panel regime where smooth quadrature
% on the fine mesh ceases to be enough.
n_far = size(targets_far,2);
err_near = err_rcip(n_far+1:end);
idx_above_1em5 = find(dvec >= 1e-5);
assert(max(err_near(idx_above_1em5)) < 1e-9, ...
    'rcip-aware should be ~1e-10 down to d=1e-5 from a corner');
% The full sweep should stay better than 1e-9 for the deepest targets.
assert(max(err_near) < 1e-9, ...
    'rcip-aware should stay <1e-9 even at sub-panel distances');
fprintf('\nTest passed.\n');

function [r,d,d2] = circle(t)
    c = cos(t(:).'); s = sin(t(:).');
    r = [c; s]; d = [-s; c]; d2 = [-c; -s];
end
