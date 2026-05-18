% helmos_demo.m
%
% End-to-end demo of helmos.m: pick a geometry, solve the Helmholtz
% Dirichlet open-arc problem with a planewave incident wave, evaluate the
% scattered field on a tensor grid, and plot amplitude + log10 error
% against a refined reference. Mirrors testhelmos's plot output.
%
% Pick `name` below to switch between geometries.

clearvars; close all; format long e;

name = 'corners1';   % 'linesegment' | 'spiral' | 'corners1' | 'corners3'
ngr  = 200;          % grid resolution per axis
ppw  = 18;
ppw_ref = ppw * 2;
nsub = 30;
theta = pi/4;
zk_scale = 10*pi;    % so that L/lambda ~ 5

build_map = struct( ...
    'linesegment', @(zk,p) geo_linesegment(zk,p), ...
    'spiral',      @(zk,p) geo_spiral(zk,p), ...
    'corners1',    @(zk,p) geo_corners(1,zk,p), ...
    'corners3',    @(zk,p) geo_corners(3,zk,p));
build = build_map.(name);

cg_probe = build(1.0, ppw);
Ltot = 0;
for ie = 1:length(cg_probe.echnks); Ltot = Ltot + sum(cg_probe.echnks(ie).wts(:)); end
zk = zk_scale / Ltot;

cg = build(zk, ppw);
cg_ref = build(zk, ppw_ref);

% Build target grid spanning the chunkgraph bounding box, padded by ~30%.
rmin = min(cg); rmax = max(cg);
cx = (rmin(1)+rmax(1))/2;  cy = (rmin(2)+rmax(2))/2;
boxw = max(rmax(1)-rmin(1), rmax(2)-rmin(2)) * 1.4;
xg = linspace(cx - boxw/2, cx + boxw/2, ngr);
yg = linspace(cy - boxw/2, cy + boxw/2, ngr);
[X,Y] = meshgrid(xg, yg);
targets = [X(:).'; Y(:).'];

ubdry_fn = @(s,b,k) helmos_planewave_rhs(s,b,k,theta);
opts = struct('nsub',nsub,'tol',1e-13,'verbose',true);

geo  = struct('cg',cg,    'targets',targets, 'ubdry_fn',ubdry_fn);
geor = struct('cg',cg_ref,'targets',targets, 'ubdry_fn',ubdry_fn);

fprintf('\n== solving %s, zk=%.3f ==\n', name, zk);
out  = helmos(geo,  'd', zk, opts);
fprintf('\n== reference solve (%s, ppw=%d) ==\n', name, ppw_ref);
outr = helmos(geor, 'd', zk, opts);

% scattered field amplitude
helmos_fieldplot(out.ufield, xg, yg, 'mode','abs', 'cg',cg, ...
    'title', sprintf('%s Dirichlet |u_{scat}| (zk=%.2f)', name, zk));

% log10 error against refined reference
helmos_fieldplot(out.ufield - outr.ufield, xg, yg, 'mode','log10err', ...
    'cg',cg, 'title', sprintf('%s log10 |u_{coarse} - u_{fine}|', name));

% Print quick summary
errfield = abs(out.ufield - outr.ufield);
fprintf('\nfield max |u-uref| = %.3e\n', max(errfield));
fprintf('field median |u-uref| = %.3e\n', median(errfield));
