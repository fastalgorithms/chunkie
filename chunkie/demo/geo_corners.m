function cg = geo_corners(ncorner, zk, ppw)
%GEO_CORNERS chain of ncorner+1 line segments connecting alternating
% heights, mirroring testhelmos's icase=2 with curvetype=1 (default).
% Vertices: ncorner+2 (2 endpoints + ncorner interior corners).
% Edges: ncorner+1 (line segments between consecutive vertices).

if nargin < 1 || isempty(ncorner), ncorner = 1; end
if nargin < 3 || isempty(ppw), ppw = 10; end
pref = []; pref.k = 16;
maxchunklen = pref.k/ppw/abs(zk)*2*pi;

% Same vertex layout as testhelmos's curvetype=1 in ossetup_corners
nsing  = ncorner + 2;
ncurve = ncorner + 1;
d = 1;
z = zeros(nsing,1);
z(1:2:nsing) = (-ncorner/2-1)*d - 0.5i + (1:2:nsing)*d;
z(2:2:nsing) = (-ncorner/2-1)*d + (2:2:nsing)*d + 0.5i;

verts = [real(z).'; imag(z).'];
edge2verts = [1:ncurve; 2:ncurve+1];

cparams = struct('ta',0,'tb',1,'maxchunklen',maxchunklen,'nover',2);
cg = chunkgraph(verts, edge2verts, [], cparams, pref);
cg = balance(cg);
end
