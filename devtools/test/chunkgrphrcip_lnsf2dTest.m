chunkgrphrcip_lnsf2dTest0();

function chunkgrphrcip_lnsf2dTest0()
%CHUNKGRPHRCIP_LNSF2DTEST
%
% Exterior Dirichlet (no-slip + isothermal) problem for the 2D linearized
% Navier-Stokes-Fourier system on a domain with corners, solved with the
% chunkgraph + RCIP machinery.
%
% RCIP requires the delta-function part of the system to be exactly the
% identity, so the combined kernel is scaled by 2:
%     (I + 2(D + eta S)) psi = 2 f,   u = (D + eta S) psi.
%
%   (1)  large-nu test (all wavenumbers O(1), boundary layer thicker than
%        the chunks): plain ggq quadrature suffices; MMS with interior
%        point sources; checks the rcip_matcorrect hook is a no-op here
%   (2)  moderate-nu test: the vertex-adjacent chunks are long compared
%        to delta, so the coarse RCIP levels need the window-quadrature
%        correction; solves with and without the correction (outer matrix
%        corrected with vertex-star pairs skipped, RCIP levels corrected
%        through the rcip_matcorrect hook)
%   (3)  realistic nu (delta ~ 3e-3, chunks 100-200 delta): the plain ggq
%        rules fail at ~1e-4; corrected solve reaches ~1e-7 on the
%        default grid
%   (3b) outer-grid self-convergence at realistic nu
%   (3c) air-like ratio (|kv|/|ka| = 1e4, delta ~ 1.4e-4): 10 digits on
%        the resolved grid
%
% see also CHNK.LNSF2D.RCIPCORRECT, CHNK.RCIP.RCOMPCHUNK

iseed = 8675309;
rng(iseed);

% ---------------- geometry: curved triangle chunkgraph ----------------
verts = [1, -0.7, -0.2; 0, 0.8, -0.9];
edge2verts = [-1, 1, 0; 0, -1, 1; 1, 0, -1];
edge2verts = sparse(edge2verts);
amp = 0.15; frq = 2.5;
fchnks = cell(1,3);
for icurve = 1:3
    fchnks{icurve} = @(t) sinearc(t, amp, frq);
end
cgrph = chunkgraph(verts, edge2verts, fchnks);
bdrypts = reshape(cgrph.r, 2, []);

% interior sources and exterior targets for MMS
srcpts = []; srcpts.r = [0.05, -0.15; -0.05, 0.10];
strengths = randn(6,1) + 1i*randn(6,1);
tt = 2*pi*(0:7)/8;
targs = [2.5*cos(tt); 2.5*sin(tt)];

% ---------------- (1) large nu ----------------
popts = []; popts.mu = 0.4; popts.muB = 0.5; popts.kap = 0.625;
prm = chnk.lnsf2d.params(popts);
eta = prm.mu;
fprintf('large nu:  ka = %.3f%+.3fi, |kh| = %.3f, |kv| = %.3f, delta = %.3f\n', ...
    real(prm.ka), imag(prm.ka), abs(prm.kh), abs(prm.kv), 1/imag(prm.kh));

fkern = kernel(@(s,t) 2*chnk.lnsf2d.kern(prm, s, t, 'c', eta));
fkern.opdims = [3 3];
fkern.sing = 'log';
fkerneval = kernel(@(s,t) chnk.lnsf2d.kern(prm, s, t, 'c', eta));
fkerneval.opdims = [3 3];
fkerneval.sing = 'log';

uex = @(pts) chnk.lnsf2d.kern(prm, srcpts, struct('r', pts), 's')*strengths;
utrue = uex(targs);

opts = []; opts.rcip = true; opts.nonsmoothonly = false;
sysmat = chunkermat(cgrph, fkern, opts);
sys = eye(size(sysmat,1)) + sysmat;
rhs = 2*uex(bdrypts);
psi = sys \ rhs;
optse = []; optse.usesmooth = false;
usol = chunkerkerneval(cgrph, fkerneval, psi, targs, optse);
relerr = max(abs(usol - utrue))/max(abs(utrue));
fprintf('large nu:  corner MMS relerr = %.3e\n', relerr);
assert(relerr < 1e-6);

% hook plumbing: at large nu the correction must be a no-op
opts2 = opts;
opts2.rcip_matcorrect = chnk.lnsf2d.rcipcorrect(prm, eta);
sysmat2 = chunkermat(cgrph, fkern, opts2);
fprintf('large nu:  hook no-op check, ||diff|| = %.3e\n', ...
    norm(sysmat2 - sysmat, 'fro'));
assert(norm(sysmat2 - sysmat, 'fro') == 0);

% ---------------- (2) moderate nu ----------------
popts = []; popts.mu = 0.004; popts.muB = 0.005; popts.kap = 0.00625;
prm2 = chnk.lnsf2d.params(popts);
eta2 = prm2.mu;
delta2 = 1/imag(prm2.kh);
fprintf('moderate nu: |kv|/|ka| = %.1f, delta = %.4f\n', ...
    abs(prm2.kv)/abs(prm2.ka), delta2);

fkern2 = kernel(@(s,t) 2*chnk.lnsf2d.kern(prm2, s, t, 'c', eta2));
fkern2.opdims = [3 3];
fkern2.sing = 'log';
fkerneval2 = kernel(@(s,t) chnk.lnsf2d.kern(prm2, s, t, 'c', eta2));
fkerneval2.opdims = [3 3];
fkerneval2.sing = 'log';
uex2 = @(pts) chnk.lnsf2d.kern(prm2, srcpts, struct('r', pts), 's')*strengths;
utrue2 = uex2(targs);
rhs2 = 2*uex2(bdrypts);

merged = merge(cgrph.echnks);
lens = chunklen(merged);
flags = chnk.lnsf2d.flagchunks(merged, prm2, 2.0);
fprintf('moderate nu: max chunklen/delta = %.2f, flagged chunks: %d of %d\n', ...
    max(lens)/delta2, numel(flags), merged.nch);

res2 = zeros(1,2);
for iv = 1:2
    optsv = []; optsv.rcip = true; optsv.nonsmoothonly = false;
    if iv == 2
        copts = []; copts.tau = 2.0;
        optsv.rcip_matcorrect = chnk.lnsf2d.rcipcorrect(prm2, eta2, copts);
    end
    [sysmat, ~, rcipsav] = chunkermat(cgrph, fkern2, optsv);
    if iv == 2
        % outer correction: skip pairs owned by the rcip vertex blocks
        skipsets = cell(numel(rcipsav), 1);
        for v = 1:numel(rcipsav)
            if ~isempty(rcipsav{v}) && isfield(rcipsav{v}, 'starind')
                skipsets{v} = unique(ceil(rcipsav{v}.starind/3));
            end
        end
        copts2 = copts; copts2.scale = 2.0; copts2.skipsets = skipsets;
        sysmat = chnk.lnsf2d.nearcorrect(sysmat, merged, prm2, eta2, copts2);
    end
    sys = eye(size(sysmat,1)) + sysmat;
    psi = sys \ rhs2;
    usol = chunkerkerneval(cgrph, fkerneval2, psi, targs, optse);
    res2(iv) = max(abs(usol - utrue2))/max(abs(utrue2));
end
fprintf('moderate nu: plain rcip = %.3e   corrected rcip = %.3e\n', ...
    res2(1), res2(2));
assert(res2(2) < 1e-5);
assert(res2(2) < 2*res2(1) + 1e-9);

% ---------------- (3) realistic nu: the strain test ----------------
% delta ~ 3e-3, chunks are 100-200 delta long: the plain ggq self/near
% rules fail at ~1e-4, and the top ~7 levels of each rcip recursion are
% flagged. The outer matrix is corrected with the exact rcip starind
% blocks skipped (they are owned by inv(R) - I).
popts = []; popts.mu = 4e-6; popts.muB = 5e-6; popts.kap = 6.25e-6;
prm3 = chnk.lnsf2d.params(popts);
eta3 = prm3.mu;
delta3 = 1/imag(prm3.kh);
fprintf('realistic nu: |kv|/|ka| = %.0f, delta = %.2e, maxlen/delta = %.0f\n', ...
    abs(prm3.kv)/abs(prm3.ka), delta3, max(lens)/delta3);

fkern3 = kernel(@(s,t) 2*chnk.lnsf2d.kern(prm3, s, t, 'c', eta3));
fkern3.opdims = [3 3];
fkern3.sing = 'log';
fkerneval3 = kernel(@(s,t) chnk.lnsf2d.kern(prm3, s, t, 'c', eta3));
fkerneval3.opdims = [3 3];
fkerneval3.sing = 'log';
uex3 = @(pts) chnk.lnsf2d.kern(prm3, srcpts, struct('r', pts), 's')*strengths;
utrue3 = uex3(targs);
rhs3 = 2*uex3(bdrypts);

res3 = zeros(1,2);
for iv = 1:2
    optsv = []; optsv.rcip = true; optsv.nonsmoothonly = false;
    if iv == 2
        copts = []; copts.tau = 2.0;
        optsv.rcip_matcorrect = chnk.lnsf2d.rcipcorrect(prm3, eta3, copts);
    end
    [sysmat, ~, rcipsav] = chunkermat(cgrph, fkern3, optsv);
    if iv == 2
        % exact skip sets from the rcip starind blocks (ndim = 3)
        skipsets = cell(numel(rcipsav), 1);
        for v = 1:numel(rcipsav)
            if ~isempty(rcipsav{v}) && isfield(rcipsav{v}, 'starind')
                skipsets{v} = unique(ceil(rcipsav{v}.starind/3));
            end
        end
        copts2 = copts; copts2.scale = 2.0; copts2.skipsets = skipsets;
        sysmat = chnk.lnsf2d.nearcorrect(sysmat, merged, prm3, eta3, copts2);
    end
    sys = eye(size(sysmat,1)) + sysmat;
    psi = sys \ rhs3;
    usol = chunkerkerneval(cgrph, fkerneval3, psi, targs, optse);
    res3(iv) = max(abs(usol - utrue3))/max(abs(utrue3));
end
fprintf('realistic nu: plain rcip = %.3e   corrected rcip = %.3e\n', ...
    res3(1), res3(2));
assert(res3(2) < 1e-6);
assert(res3(1)/res3(2) > 1e2);   % the correction should win big here

% ---------------- (3b) outer-grid self-convergence ---------------------
% refine the outer grid at fixed physics: the corrected corner solve
% converges like a smooth-boundary discretization even though the chunks
% remain ~80 delta long, i.e. the grid requirement is set by geometry and
% data, not by delta
for cscale = [0.5, 0.25]
    cparams = cell(1,3);
    for icurve = 1:3
        cparams{icurve} = []; cparams{icurve}.maxchunklen = cscale;
    end
    cgrph2 = chunkgraph(verts, edge2verts, fchnks, cparams);
    merged2 = merge(cgrph2.echnks);
    bdrypts2 = reshape(cgrph2.r, 2, []);
    optsv = []; optsv.rcip = true; optsv.nonsmoothonly = false;
    copts = []; copts.tau = 2.0;
    optsv.rcip_matcorrect = chnk.lnsf2d.rcipcorrect(prm3, eta3, copts);
    [sysmat, ~, rcipsav] = chunkermat(cgrph2, fkern3, optsv);
    skipsets = cell(numel(rcipsav), 1);
    for v = 1:numel(rcipsav)
        if ~isempty(rcipsav{v}) && isfield(rcipsav{v}, 'starind')
            skipsets{v} = unique(ceil(rcipsav{v}.starind/3));
        end
    end
    copts2 = copts; copts2.scale = 2.0; copts2.skipsets = skipsets;
    sysmat = chnk.lnsf2d.nearcorrect(sysmat, merged2, prm3, eta3, copts2);
    sys = eye(size(sysmat,1)) + sysmat;
    psi2 = sys \ (2*uex3(bdrypts2));
    usol = chunkerkerneval(cgrph2, fkerneval3, psi2, targs, optse);
    errc = max(abs(usol - utrue3))/max(abs(utrue3));
    fprintf('realistic nu: maxchunklen = %.2f (nch = %d): corrected = %.3e\n', ...
        cscale, merged2.nch, errc);
end
% on the resolved grid the corrected corner solve reaches the smooth-
% boundary floor even though the chunks are still ~80 delta long
assert(errc < 1e-9);

% ---------------- (3c) air-like ratio at fixed grids -------------------
% |kv|/|ka| = 1e4, delta ~ 1.4e-4. If the grid requirements are set by
% geometry/data only (delta-independent), the fine grid (nch = 48, chunks
% still ~900 delta) should hold its floor and the coarse grid should
% reproduce its ~7e-7 unchanged. This also stress-tests the integration
% of the "almost-pv" static tails, whose graded-piece count grows only
% like log(L/delta).
%
% Grid note: the outer-grid discretization error (density interpolation,
% independent of delta and untouchable by quadrature or rcip settings)
% is 4.5e-7 / 6.9e-10 / 6.8e-11 at maxchunklen 0.5 / 0.25 / 0.125 for
% this geometry, so 10 digits requires the 0.125 grid.
popts = []; popts.mu = 1e-8; popts.muB = 1.25e-8; popts.kap = 1.5625e-8;
prm4 = chnk.lnsf2d.params(popts);
eta4 = prm4.mu;
% thermal-density conjugation: balances the eta_T/kappa ~ delta^{-2}
% block imbalance so that all inversions (backslash and the rcip
% recursion) are well conditioned uniformly in delta
tsc4 = abs(prm4.kh);
fprintf('air-like: |kv|/|ka| = %.0f, delta = %.2e\n', ...
    abs(prm4.kv)/abs(prm4.ka), 1/imag(prm4.kh));
fkern4 = kernel(@(s,t) 2*chnk.lnsf2d.kern(prm4, s, t, 'c', eta4, tsc4));
fkern4.opdims = [3 3];
fkern4.sing = 'log';
fkerneval4 = kernel(@(s,t) chnk.lnsf2d.kern(prm4, s, t, 'c', eta4));
fkerneval4.opdims = [3 3];
fkerneval4.sing = 'log';
uex4 = @(pts) chnk.lnsf2d.kern(prm4, srcpts, struct('r', pts), 's')*strengths;
utrue4 = uex4(targs);

cparams = cell(1,3);
for icurve = 1:3
    cparams{icurve} = []; cparams{icurve}.maxchunklen = 0.125;
end
cgrph4 = chunkgraph(verts, edge2verts, fchnks, cparams);

grids = {cgrph, cgrph4};    % coarse (nch 12) and fine (nch 48)
glabels = {'coarse', 'fine  '};
err4 = zeros(1,2);
for ig = 1:2
    cg = grids{ig};
    mg = merge(cg.echnks);
    optsv = []; optsv.rcip = true; optsv.nonsmoothonly = false;
    copts = []; copts.tau = 2.0; copts.tscale = tsc4;
    optsv.rcip_matcorrect = chnk.lnsf2d.rcipcorrect(prm4, eta4, copts);
    [sysmat, ~, rcipsav] = chunkermat(cg, fkern4, optsv);
    skipsets = cell(numel(rcipsav), 1);
    for v = 1:numel(rcipsav)
        if ~isempty(rcipsav{v}) && isfield(rcipsav{v}, 'starind')
            skipsets{v} = unique(ceil(rcipsav{v}.starind/3));
        end
    end
    copts2 = copts; copts2.scale = 2.0; copts2.skipsets = skipsets;
    copts2.tscale = tsc4;
    sysmat = chnk.lnsf2d.nearcorrect(sysmat, mg, prm4, eta4, copts2);
    sys = eye(size(sysmat,1)) + sysmat;
    rhs4 = 2*uex4(reshape(cg.r, 2, []));
    rhs4(3:3:end) = rhs4(3:3:end)/tsc4;       % conjugated rhs
    psi4 = sys \ rhs4;
    psi4(3:3:end) = psi4(3:3:end)*tsc4;       % physical density
    usol = chunkerkerneval(cg, fkerneval4, psi4, targs, optse);
    err4(ig) = max(abs(usol - utrue4))/max(abs(utrue4));
    fprintf('air-like: %s grid corrected = %.3e\n', glabels{ig}, err4(ig));
end
% the corrected quadrature + rcip floor is below 1e-10 at air (measured
% 6.8e-11 on this grid); the coarse-grid result is pure outer
% discretization error
assert(err4(2) < 1e-10);

end

% ------------------------------------------------------------------
function [r,d,d2] = sinearc(t,amp,frq)
xs = t;
ys = amp*sin(frq*t);
xp = ones(size(t));
yp = amp*frq*cos(frq*t);
xpp = zeros(size(t));
ypp = -frq*frq*amp*sin(frq*t);
r = [(xs(:)).'; (ys(:)).'];
d = [(xp(:)).'; (yp(:)).'];
d2 = [(xpp(:)).'; (ypp(:)).'];
end
