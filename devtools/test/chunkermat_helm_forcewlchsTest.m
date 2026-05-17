function chunkermat_helm_forcewlchsTest()
%CHUNKERMAT_HELM_FORCEWLCHSTEST  Verify chunkermat opts.forcewlchs dispatches
%   correctly to chnk.kernsplit Helmholtz self/adj corrections (s, d, sp).
%
%   On a smooth closed curve, the wLCHS quadrature and chunkermat's default
%   GGQ quadrature should agree in their action on a smooth density.
%   Comparison tolerance ~1e-12 (limited by GGQ).

zk = 1.1;
cparams = []; cparams.ta = 0; cparams.tb = 2*pi; cparams.ifclosed = true;
cparams.maxchunklen = 2*pi/8;
pref = []; pref.k = 16;
chnkr = chunkerfunc(@(t) circle_param(t), cparams, pref);
N = chnkr.npt;

% Smooth density: complex exponential along the boundary.
theta = atan2(chnkr.r(2,:), chnkr.r(1,:));
sig = exp(2i*theta(:));

types = {'s','d','sp','dp'};
% Tolerances reflect GGQ quadrature accuracy (the reference); wlchs is
% typically more accurate. dp is hypersingular so GGQ tolerance is loose.
tols  = [1e-11 1e-11 1e-10 1e-6];
for it = 1:numel(types)
    t = types{it};
    K = kernel('helm', t, zk);
    A_ref = chunkermat(chnkr, K);
    A_wl  = chunkermat(chnkr, K, struct('forcewlchs',true));
    diff = norm((A_ref - A_wl) * sig, 'inf');
    fprintf('  helm %-3s: ||(A_ggq - A_wlchs)*sig||_inf = %.3e\n', t, diff);
    assert(diff < tols(it), 'Helmholtz wlchs %s mismatch %.3e exceeds %.1e', ...
        t, diff, tols(it));
end

% Scalar-multiplied kernel: 3.7 * S
kern_scaled = 3.7 * kernel('helm','s',zk);
A_ref = chunkermat(chnkr, kern_scaled);
A_wl  = chunkermat(chnkr, kern_scaled, struct('forcewlchs',true));
diff = norm((A_ref - A_wl) * sig, 'inf');
fprintf('  3.7*s   : ||(A_ggq - A_wlchs)*sig||_inf = %.3e\n', diff);
assert(diff < 1e-11);

% Combined kernels ('c' = c1*D + c2*S, 'sc' = c1*S' + c2*S) are also
% supported in the layout dispatch, but their self block requires RCIP
% corrections for full accuracy -- they are exercised in the helmopensurface
% tests of the open-arc Helmholtz PR.

fprintf('PASS: chunkermat forcewlchs helm dispatch verified for s/d/sp/dp/scaled\n');
end

function [r, d, d2] = circle_param(t)
ct = cos(t); st = sin(t);
r  = [ct(:).';  st(:).'];
d  = [-st(:).'; ct(:).'];
d2 = [-ct(:).'; -st(:).'];
end
