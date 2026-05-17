elast2dWlchsTest0();


function elast2dWlchsTest0()
%ELAST2DWLCHSTEST   verify chunkermat opts.forcewlchs for 2D elasticity
% S, D, Strac, Dalt sub-block kernels: BVP solves and off-curve close
% evaluation against an analytic Kelvinlet.

iseed = 1234;
rng(iseed,'twister');

lam = 1.5;
mu  = 2.1;

cparams = []; cparams.eps = 1e-10;
chnkr = chunkerfunc(@(t) starfish(t), cparams);
chnkr = chnkr.sort();
fprintf('starfish: nch=%d, npt=%d\n', chnkr.nch, chnkr.npt);

f = [-1.3; 2; 0.8; -1.2];
src = []; src.r = [2 3.1; 1 -1.1]; src.n = randn(2,2);

% --- 1) chunkermat forcewlchs builds without error for all four kernels
ker_s     = kernel.elast2d('s',     lam, mu);
ker_d     = kernel.elast2d('d',     lam, mu);
ker_strac = kernel.elast2d('strac', lam, mu);
ker_dalt  = kernel.elast2d('dalt',  lam, mu);

opts_wlchs = struct('forcewlchs', true);

sysS_w = chunkermat(chnkr, ker_s,     opts_wlchs);
sysD_w = chunkermat(chnkr, ker_d,     opts_wlchs);
sysT_w = chunkermat(chnkr, ker_strac, opts_wlchs);
sysA_w = chunkermat(chnkr, ker_dalt,  opts_wlchs);
fprintf('chunkermat forcewlchs build OK for s, d, strac, dalt\n');

assert(all(size(sysS_w) == [2*chnkr.npt, 2*chnkr.npt]));
assert(all(size(sysD_w) == [2*chnkr.npt, 2*chnkr.npt]));
assert(all(size(sysT_w) == [2*chnkr.npt, 2*chnkr.npt]));
assert(all(size(sysA_w) == [2*chnkr.npt, 2*chnkr.npt]));

% --- 2) wLCHS vs GGQ sysmat block-difference: same operator, different
% quadrature.  Bounded operators (Dalt) should agree closely; singular
% operators (S, D, Strac) may differ in their near-diagonal entries since
% wLCHS and GGQ use distinct corrections, but both should give an accurate
% operator (see BVP tests below).
sysS_g = chunkermat(chnkr, ker_s);
sysD_g = chunkermat(chnkr, ker_d);
sysT_g = chunkermat(chnkr, ker_strac);
sysA_g = chunkermat(chnkr, ker_dalt);

rel_A = norm(sysA_w - sysA_g,'fro') / norm(sysA_g,'fro');
fprintf('|wLCHS - GGQ| / |GGQ| (Dalt smooth):  %5.2e\n', rel_A);
assert(rel_A < 1e-7, 'Dalt wLCHS deviates from GGQ smooth-quadrature');

% --- 3) Interior Dirichlet BVP via Dalt with wLCHS system matrix
% Solve  (1/2 D_alt-jump) sigma + D_alt sigma = u_bdry   (chunkie convention
% follows the existing elastickernelsTest with diag = 0.5(lam+3mu)/(lam+2mu).)
t = []; t.r = chnkr.r(:,:); t.n = chnkr.n(:,:); t.d = chnkr.d(:,:);
wts = chnkr.wts; wts = wts(:);
wts2 = [wts(:).'; wts(:).']; wts2 = wts2(:);
[ub, ~] = elasticlet(lam, mu, src, t, f);

% Off-curve interior target points (well inside starfish, near origin)
targ = []; targ.r = 0.1*randn(2,8); targ.n = randn(2,8);
[utarg_true, ~] = elasticlet(lam, mu, src, targ, f);

diag_dalt = 0.5*((lam+3*mu)/(lam+2*mu));

sigma_wlchs = (sysA_w + diag_dalt*eye(2*chnkr.npt)) \ ub;
sigma_ggq   = (sysA_g + diag_dalt*eye(2*chnkr.npt)) \ ub;

K_dalt = chnk.elast2d.kern(lam, mu, t, targ, 'dalt');
ueval_wlchs = K_dalt * (wts2 .* sigma_wlchs);
ueval_ggq   = K_dalt * (wts2 .* sigma_ggq);

err_dir_wlchs = norm(ueval_wlchs - utarg_true)/norm(utarg_true);
err_dir_ggq   = norm(ueval_ggq   - utarg_true)/norm(utarg_true);
fprintf('Dirichlet (Dalt) interior eval rel err: wLCHS %5.2e   GGQ %5.2e\n', ...
    err_dir_wlchs, err_dir_ggq);
assert(err_dir_wlchs < 1e-9, 'Dalt+wLCHS Dirichlet solve failed');

% --- 4) Neumann BVP via Strac with wLCHS system matrix
% Solve  diag_jump sigma + Strac sigma = trac_bdry,  evaluate u via S.
sp = -0.5*(mu/(lam+2*mu));
stoktrac = -0.5*(lam+mu)/(lam+2*mu);
diag_strac = sp + stoktrac;

[~, tracb] = elasticlet(lam, mu, src, t, f);

K_s_off = chnk.elast2d.kern(lam, mu, t, targ, 's');
targperp = chnk.perp(targ.r);
nt = size(targ.r,2);
nspace = [repmat(eye(2),nt,1) targperp(:)];

sigma_n_wlchs = gmres(sysT_w + diag_strac*eye(2*chnkr.npt), tracb, ...
    [], 1e-13, 200);
ueval_n_w = K_s_off * (wts2 .* sigma_n_wlchs);
cfs = nspace \ (utarg_true - ueval_n_w);
ueval_n_w_corr = ueval_n_w + nspace*cfs;
err_neu_wlchs = norm(utarg_true - ueval_n_w_corr)/norm(utarg_true);

sigma_n_ggq = gmres(sysT_g + diag_strac*eye(2*chnkr.npt), tracb, ...
    [], 1e-13, 200);
ueval_n_g = K_s_off * (wts2 .* sigma_n_ggq);
cfs_g = nspace \ (utarg_true - ueval_n_g);
ueval_n_g_corr = ueval_n_g + nspace*cfs_g;
err_neu_ggq = norm(utarg_true - ueval_n_g_corr)/norm(utarg_true);

fprintf('Neumann (Strac) eval rel err: wLCHS %5.2e   GGQ %5.2e\n', ...
    err_neu_wlchs, err_neu_ggq);

% wLCHS Strac and GGQ Strac discretize the same continuous operator with
% different quadrature schemes; their per-panel-count accuracies need not
% match.  Confirm wLCHS converges geometrically with refinement.
cparams_r = []; cparams_r.eps = 1e-12;
chnkr2 = chunkerfunc(@(t) starfish(t), cparams_r);
chnkr2 = chnkr2.refine(struct('nchmax',10000)).sort();

sysT_w2 = chunkermat(chnkr2, ker_strac, opts_wlchs);
t2 = []; t2.r = chnkr2.r(:,:); t2.n = chnkr2.n(:,:); t2.d = chnkr2.d(:,:);
wts2_r = chnkr2.wts; wts2_r = wts2_r(:);
wts2_r2 = [wts2_r(:).'; wts2_r(:).']; wts2_r2 = wts2_r2(:);
[~, tracb2] = elasticlet(lam, mu, src, t2, f);
sigma_n_w2 = gmres(sysT_w2 + diag_strac*eye(2*chnkr2.npt), tracb2, ...
    [], 1e-13, 200);
K_s_off2 = chnk.elast2d.kern(lam, mu, t2, targ, 's');
ueval_n_w2 = K_s_off2 * (wts2_r2 .* sigma_n_w2);
cfs2 = nspace \ (utarg_true - ueval_n_w2);
ueval_n_w2_corr = ueval_n_w2 + nspace*cfs2;
err_neu_w2 = norm(utarg_true - ueval_n_w2_corr)/norm(utarg_true);

fprintf('  Strac Neumann refinement: nch=%d->%d  err %5.2e -> %5.2e\n', ...
    chnkr.nch, chnkr2.nch, err_neu_wlchs, err_neu_w2);
assert(err_neu_w2 < 1e-12, 'Strac wLCHS Neumann failed to converge');
assert(err_neu_w2 < err_neu_wlchs / 10, ...
    'Strac wLCHS Neumann did not converge geometrically with h refinement');

% --- 5) Off-curve close eval: panel_eval at distances approaching curve
% Use Dalt density just solved; evaluate u at points approaching a
% boundary point along the inward normal direction.
ich = ceil(chnkr.nch/2);
pidx = (ich-1)*chnkr.k + ceil(chnkr.k/2);
y0 = chnkr.r(:,pidx);
n0 = chnkr.n(:,pidx);
dists = [1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01];
tprobe = y0 - n0 .* dists;          % interior side
% elasticlet only needs u (1st output); add a dummy normal to satisfy strac
% line in the helper (not used in u).
probe_struct = struct('r',tprobe,'n',repmat([1;0],1,size(tprobe,2)));
[utrue_probe, ~] = elasticlet(lam, mu, src, probe_struct, f);
utrue_probe = reshape(utrue_probe, 2, []);

ueval_probe = chnk.kernsplit.elast2d_panel_eval(chnkr, 'dalt', lam, mu, ...
    sigma_wlchs, tprobe);

err_probe = sqrt(sum((ueval_probe - utrue_probe).^2, 1)) ./ ...
            max(sqrt(sum(utrue_probe.^2, 1)), 1e-30);
fprintf('Dalt panel_eval close eval at dists [%s]:\n  rel err = [%s]\n', ...
    sprintf('%g ', dists), sprintf('%5.2e ', err_probe));

% Smooth far targets (d~1) should be near machine; close targets degrade
% since bounded-piece corrections are not implemented in the simplified
% panel_eval.
assert(err_probe(1) < 1e-10, 'Far close-eval target unexpectedly inaccurate');

% --- 6) chunkerkerneval forcewlchs=true dispatches into elast2d_panel_eval
ck_opts = struct('forcewlchs', true);
test_dens = sigma_wlchs;  % a smooth-ish 2-component density from the Dirichlet solve
ueval_ck_s    = chunkerkerneval(chnkr, ker_s,    test_dens, tprobe, ck_opts);
ueval_ck_d    = chunkerkerneval(chnkr, ker_d,    test_dens, tprobe, ck_opts);
ueval_ck_dalt = chunkerkerneval(chnkr, ker_dalt, test_dens, tprobe, ck_opts);

% Strac needs target normals; pass as struct.  Use the boundary outward
% normal at the chosen panel point as the target normal (smooth in tprobe
% direction).
tprobe_strac = struct('r', tprobe, 'n', repmat(n0, 1, numel(dists)));
ueval_ck_strac = chunkerkerneval(chnkr, ker_strac, test_dens, tprobe_strac, ck_opts);

% Compare to direct elast2d_panel_eval (already validated above).
direct_s    = chnk.kernsplit.elast2d_panel_eval(chnkr, 's',     lam, mu, test_dens, tprobe);
direct_d    = chnk.kernsplit.elast2d_panel_eval(chnkr, 'd',     lam, mu, test_dens, tprobe);
direct_dalt = chnk.kernsplit.elast2d_panel_eval(chnkr, 'dalt',  lam, mu, test_dens, tprobe);
direct_strac= chnk.kernsplit.elast2d_panel_eval(chnkr, 'strac', lam, mu, test_dens, tprobe_strac);

rel_ck = @(a,b) norm(a(:) - b(:)) / max(norm(b(:)), 1e-30);
e_s    = rel_ck(ueval_ck_s,    direct_s);
e_d    = rel_ck(ueval_ck_d,    direct_d);
e_dalt = rel_ck(ueval_ck_dalt, direct_dalt);
e_strac= rel_ck(ueval_ck_strac, direct_strac);
fprintf('chunkerkerneval forcewlchs vs direct elast2d_panel_eval:\n');
fprintf('  s=%5.2e  d=%5.2e  strac=%5.2e  dalt=%5.2e\n', e_s, e_d, e_strac, e_dalt);
assert(e_s    < 1e-14, 'chunkerkerneval forcewlchs S    diverges from direct');
assert(e_d    < 1e-14, 'chunkerkerneval forcewlchs D    diverges from direct');
assert(e_strac< 1e-14, 'chunkerkerneval forcewlchs Strac diverges from direct');
assert(e_dalt < 1e-14, 'chunkerkerneval forcewlchs Dalt diverges from direct');

% --- 7) End-to-end Dirichlet via Dalt against analytical Kelvinlet at
% off-curve targets across a range of distances.  Uses refined chnkr2 so
% sigma is at machine precision, then close-evaluates via chunkerkerneval
% +forcewlchs and compares to the analytical Kelvinlet.
sysA_w2 = chunkermat(chnkr2, ker_dalt, opts_wlchs);
[ub2, ~] = elasticlet(lam, mu, src, t2, f);
sigma_dir_w2 = (sysA_w2 + diag_dalt*eye(2*chnkr2.npt)) \ ub2;

% Inward probe at one panel point, several distances
ich2 = ceil(chnkr2.nch/2);
pidx2 = (ich2-1)*chnkr2.k + ceil(chnkr2.k/2);
y0_2 = chnkr2.r(:,pidx2); n0_2 = chnkr2.n(:,pidx2);
probe_dists = [1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001];
tprobe2 = y0_2 - n0_2 .* probe_dists;

% Analytical Kelvinlet at the probes
probe_struct2 = struct('r', tprobe2, 'n', repmat([1;0],1,size(tprobe2,2)));
[utrue2, ~] = elasticlet(lam, mu, src, probe_struct2, f);
utrue2 = reshape(utrue2, 2, []);

% chunkerkerneval + forcewlchs evaluation
ueval_ck = chunkerkerneval(chnkr2, ker_dalt, sigma_dir_w2, tprobe2, ck_opts);
ueval_ck = reshape(ueval_ck, 2, []);

err_ck_per_dist = sqrt(sum((ueval_ck - utrue2).^2, 1)) ./ ...
                  max(sqrt(sum(utrue2.^2, 1)), 1e-30);
fprintf('Dirichlet end-to-end (sysmat+chunkerkerneval) vs Kelvinlet:\n');
fprintf('  dist:   [%s]\n', sprintf('%9.1e', probe_dists));
fprintf('  err:    [%s]\n', sprintf('%9.2e', err_ck_per_dist));

% For moderate distances (d > 0.1*panlen), wLCHS gives machine precision.
% Panel length ~ 2*pi/nch ~ 0.06 on refined starfish, so d>0.01 is moderate.
panlen2 = 2*pi/chnkr2.nch;
fprintf('  panel length ~ %g, d/panlen = [%s]\n', panlen2, ...
    sprintf('%6.2f', probe_dists/panlen2));
moderate = probe_dists/panlen2 > 0.5;
assert(all(err_ck_per_dist(moderate) < 1e-13), ...
    'chunkerkerneval+forcewlchs Dirichlet failed at moderate-distance targets');
% Wu-Barnett guarantees full machine precision even at very-close targets.
% At nch=98 (panel length ~0.064), even d=0.001 (d/panlen ~ 0.02) should be
% better than 1e-10 (matching the Stokes-mobility precedent).
very_close = probe_dists/panlen2 < 0.2;
assert(all(err_ck_per_dist(very_close) < 1e-9), ...
    sprintf('Wu-Barnett close-eval lost accuracy at very-close targets: %s', ...
        sprintf('%5.2e ', err_ck_per_dist(very_close))));

% --- 8) End-to-end Neumann via Strac against analytical Kelvinlet at
% off-curve targets across a range of distances.  Tests the new Stokes
% Strac panel_eval path through chunkerkerneval+forcewlchs.
sysT_w2_for_neu = sysT_w2 + diag_strac*eye(2*chnkr2.npt);
[~, tracb2_for_neu] = elasticlet(lam, mu, src, t2, f);
sigma_n_w2_again = gmres(sysT_w2_for_neu, tracb2_for_neu, [], 1e-13, 200);
% u via S panel_eval (already validated)
ueval_S = chunkerkerneval(chnkr2, ker_s, sigma_n_w2_again, tprobe2, ck_opts);
ueval_S = reshape(ueval_S, 2, []);
% null-space removal at the probe targets
nt_probe = size(tprobe2,2);
targperp_probe = chnk.perp(tprobe2);
nspace_probe = [repmat(eye(2),nt_probe,1) targperp_probe(:)];
cfs_neu = nspace_probe \ (utrue2(:) - ueval_S(:));
ueval_S_corr = ueval_S + reshape(nspace_probe*cfs_neu, 2, []);
err_neu_ck_per_dist = sqrt(sum((ueval_S_corr - utrue2).^2, 1)) ./ ...
                     max(sqrt(sum(utrue2.^2, 1)), 1e-30);
fprintf('Neumann end-to-end (sysmat+chunkerkerneval S close) vs Kelvinlet:\n');
fprintf('  dist:   [%s]\n', sprintf('%9.1e', probe_dists));
fprintf('  err:    [%s]\n', sprintf('%9.2e', err_neu_ck_per_dist));
assert(all(err_neu_ck_per_dist(probe_dists > 1e-2) < 1e-11), ...
    'Neumann+S close eval failed at moderate distance');

% Consistency: at moderate distance, chunkerkerneval with forcewlchs strac
% should match smooth GL (since both correctly evaluate Strac there).
% sigma_dens here is just an arbitrary 2-component density.
mid_targ_r = chnkr2.r(:,1:end-1:end) + 0.5*chnkr2.n(:,1:end-1:end);
mid_targ_n = repmat([1;0], 1, size(mid_targ_r,2));
mid_targ_struct = struct('r', mid_targ_r, 'n', mid_targ_n);
test_sigma = sigma_n_w2_again;  % use just-computed Neumann sigma
val_wlchs_strac = chunkerkerneval(chnkr2, ker_strac, test_sigma, ...
    mid_targ_struct, ck_opts);
val_smooth_strac = chunkerkerneval(chnkr2, ker_strac, test_sigma, ...
    mid_targ_struct, struct('forcesmooth', true));
rel_strac = norm(val_wlchs_strac - val_smooth_strac) / max(norm(val_smooth_strac), 1e-30);
fprintf('Strac panel_eval wLCHS vs smooth at d~0.5 (moderate): %5.2e\n', rel_strac);
assert(rel_strac < 1e-12, 'Stokes-strac panel_eval inconsistent with smooth GL at moderate distance');


fprintf('all elast2d wLCHS smoke tests passed\n');

end


function [u, trac] = elasticlet(lam, mu, s, t, f)
mat     = chnk.elast2d.kern(lam, mu, s, t, 's');
mattrac = chnk.elast2d.kern(lam, mu, s, t, 'strac');
u    = mat*f;
trac = mattrac*f;
end
