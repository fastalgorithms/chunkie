function cg = geo_linesegment(zk, ppw, klam)
%GEO_LINESEGMENT chunkgraph of a line segment from (-1,-0.2) to (1,-0.2)
% (matching testhelmos's icase=0 placement).
%
% Inputs:
%   zk   - wavenumber
%   ppw  - points per wavelength (default 10)
%   klam - desired chunk length scaling factor; if empty, default = 1.

if nargin < 2 || isempty(ppw), ppw = 10; end
if nargin < 3 || isempty(klam), klam = 1; end

pref = []; pref.k = 16;
maxchunklen = klam * pref.k/ppw/abs(zk)*2*pi;
verts = [-1, 1; -0.2, -0.2];
edge2verts = [1; 2];
cparams = struct('ta',0,'tb',1,'maxchunklen',maxchunklen,'nover',2);
cg = chunkgraph(verts, edge2verts, [], cparams, pref);
cg = balance(cg);
end
