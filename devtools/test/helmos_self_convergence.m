% helmos_self_convergence.m
%
% Compare helmos solutions at two resolutions for both RHS types and
% several target locations. Goal: isolate whether the demo's ~1e-7
% error floor is due to RHS, target location, or just noisy convergence.

clearvars; close all; format long e;

% Pick up helmos.m + helmos_planewave_rhs + geo_*.m from chunkie/demo/.
this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, '..', '..', 'chunkie', 'demo'));

zk = 3.0;
ppw_low = 10;  ppw_hi = 30;
nsub_low = 30; nsub_hi = 30;  % keep nsub fixed; deeper RCIP loses precision

% Targets: Helsing's far one + a near-vertex sweep
targets = [ 0.17, -0.99, -0.999, -0.9999;
            0.62, -0.05, -0.005, -0.0005 ];

% --- Polynomial RHS (matches dirichletTest)
ubdry_poly = @(s,b,k) (4*s.r(1,:).^3 + 2*s.r(1,:).^2 - 3*s.r(1,:) - 1).';

% --- Planewave RHS, theta = pi/4
theta = pi/4;
ubdry_pw = @(s,b,k) helmos_planewave_rhs(s,b,k,theta);

cases = { {'poly', ubdry_poly}, {'planewave', ubdry_pw} };

for ic = 1:length(cases)
    name = cases{ic}{1}; ubdry = cases{ic}{2};

    cg_lo = geo_linesegment(zk, ppw_low);
    cg_hi = geo_linesegment(zk, ppw_hi);
    fprintf('\n[%s] npts low=%d, hi=%d\n', name, cg_lo.npt, cg_hi.npt);

    geo_lo = struct('cg',cg_lo,'targets',targets,'ubdry_fn',ubdry);
    geo_hi = struct('cg',cg_hi,'targets',targets,'ubdry_fn',ubdry);

    out_lo = helmos(geo_lo, 'd', zk, struct('nsub',nsub_low,'tol',1e-12,'verbose',false));
    out_hi = helmos(geo_hi, 'd', zk, struct('nsub',nsub_hi,'tol',1e-13,'verbose',false));

    err = abs(out_lo.ufield - out_hi.ufield);
    fprintf('  %-22s %-14s %-14s\n', 'target', '|u_lo|', 'lo-hi diff');
    for i = 1:size(targets,2)
        fprintf('  (%+.4f,%+.4f)        %.4e      %.3e\n', ...
            targets(1,i), targets(2,i), abs(out_lo.ufield(i)), err(i));
    end
end
fprintf('\nDone.\n');
