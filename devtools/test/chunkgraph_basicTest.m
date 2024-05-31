%chunkgraph_basicTest
%
% test of basic chunkgraph class constructors, methods, etc.
%

clearvars
addpaths_loc();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   .  .  .  builds a simple pentagonal chunkergraph 


verts = exp(1i*2*pi*(0:4)/5);
verts = [real(verts);imag(verts)];

endverts = [1:5; [2:5 1]];

edge2verts = [-1, 1, 0, 0, 0; ...
               0,-1, 1, 0, 0; ...
               0, 0,-1, 1, 0; ...
               0, 0, 0,-1, 1; ...
               1, 0, 0, 0,-1];

edge2verts = sparse(edge2verts);
amp = 0.5;
frq = 6;
fchnks    = {};
for icurve = 1:size(edge2verts,1)
    fchnks{icurve} = @(t) sinearc(t,amp,frq);
end

prefs = [];
[cgrph1] = chunkgraph(verts,edge2verts,fchnks,prefs);
[cgrph2] = chunkgraph(verts,endverts,fchnks,prefs);

assert(nnz(cgrph1.v2emat-cgrph2.v2emat) == 0);

cgrph1 = balance(cgrph1);

% adjacent triangles

verts = [1 0 -1 0; 0 1 0 -1]; edgesendverts = [1:3, 3, 4; 2:3, 1, 4, 1];
cg = chunkgraph(verts,edgesendverts);

% a graph with two edges in the old format

verts = [1 0 1; -1 0 1]; edge2verts = [-1 1 0; 0 -1 1];
cg1 = chunkgraph(verts,edge2verts);

edgesendverts = [1:2; 2:3];
cg2 = chunkgraph(verts,edgesendverts);

assert(nnz(cg1.v2emat-cg2.v2emat) == 0);


function [r,d,d2] = sinearc(t,amp,frq)

xs = t;
ys = amp*sin(frq*t);
xp = ones(size(t));
yp = amp*frq*cos(frq*t);
xpp = zeros(size(t));
ypp = -frq*frq*amp*sin(t);

r = [(xs(:)).'; (ys(:)).'];
d = [(xp(:)).'; (yp(:)).'];
d2 = [(xpp(:)).'; (ypp(:)).'];
end
