% helmos_smoke.m
%
% Smoke test for the helmos.m driver across geometries: line segment,
% spiral, multi-corner chain. Each case solves the Helmholtz Dirichlet
% open-arc problem with a planewave RHS at zk=10*pi/Ltot (low frequency
% to keep runtime small), and validates self-convergence by comparing
% to a 3x mesh-refined solution at fixed nsub. Far-target diff < 1e-10
% expected at low frequency.

clearvars; close all; format long e;

% Pick up helmos.m + helmos_planewave_rhs + geo_*.m from chunkie/demo/.
this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, '..', '..', 'chunkie', 'demo'));

cases = struct();
cases(1).name = 'linesegment'; cases(1).build = @(zk,ppw) geo_linesegment(zk,ppw);
cases(2).name = 'spiral';      cases(2).build = @(zk,ppw) geo_spiral(zk,ppw);
cases(3).name = 'corners1';    cases(3).build = @(zk,ppw) geo_corners(1,zk,ppw);
cases(4).name = 'corners3';    cases(4).build = @(zk,ppw) geo_corners(3,zk,ppw);

ppw_lo = 12;
ppw_hi = 36;
nsub   = 30;
theta  = pi/4;
ubdry_fn = @(s,b,k) helmos_planewave_rhs(s,b,k,theta);

opts_lo = struct('nsub',nsub,'tol',1e-12,'verbose',false);
opts_hi = struct('nsub',nsub,'tol',1e-13,'verbose',false);

for ic = 1:length(cases)
    name = cases(ic).name;
    build = cases(ic).build;

    % rough length-scale to set a wavenumber
    cg_probe = build(1.0, ppw_lo);
    % use total arc length to set zk so that L/lambda ~ 5
    Ltot = sum(cg_probe.echnks(1).wts(:));
    for ie = 2:length(cg_probe.echnks); Ltot = Ltot + sum(cg_probe.echnks(ie).wts(:)); end
    zk = 10*pi/Ltot;  % L/lambda ~ 5

    cg_lo = build(zk, ppw_lo);
    cg_hi = build(zk, ppw_hi);

    % targets: a bunch of off-curve points well away from any vertex,
    % plus a sweep approaching one vertex (the first one).
    rmin = min(cg_lo); rmax = max(cg_lo);
    cx = (rmin(1)+rmax(1))/2; cy = (rmin(2)+rmax(2))/2;
    box = max(rmax-rmin) * 1.5;
    far = [cx+box*[0.4 -0.3 0.2]; cy+box*[0.3 0.4 -0.5]];
    % near-vertex sweep: approach the first vertex perpendicular to the
    % normal of its first incident edge, to avoid landing on the curve.
    v1 = cg_lo.verts(:,1);
    e1 = cg_lo.echnks(1);
    n1 = e1.n(:,1,1); n1 = n1/norm(n1);  % normal at first node of edge 1
    dvec = 10.^(-(1:4));
    near = v1 + n1*dvec;
    targets = [far, near];

    geo_lo = struct('cg',cg_lo,'targets',targets,'ubdry_fn',ubdry_fn);
    geo_hi = struct('cg',cg_hi,'targets',targets,'ubdry_fn',ubdry_fn);

    fprintf('\n[%s] zk=%.3f, npts lo=%d hi=%d\n', name, zk, cg_lo.npt, cg_hi.npt);
    out_lo = helmos(geo_lo,'d',zk,opts_lo);
    out_hi = helmos(geo_hi,'d',zk,opts_hi);
    err = abs(out_lo.ufield - out_hi.ufield);

    fprintf('  far max diff:  %.3e\n', max(err(1:size(far,2))));
    fprintf('  near max diff: %.3e (over d in %s)\n', ...
        max(err(size(far,2)+1:end)), mat2str(dvec));

    assert(max(err(1:size(far,2))) < 1e-9, ...
        sprintf('[%s] far-target convergence < 1e-9 failed', name));
    assert(max(err) < 1e-8, ...
        sprintf('[%s] all-target convergence < 1e-8 failed', name));
end
fprintf('\nAll smoke tests passed.\n');
