tochunkgraphTest0();

function tochunkgraphTest0()
%TOCHUNKGRAPHTEST tests the routine for converting chunkers to chunkgraphs

seed = 8675309;
rng(seed);

% geometry parameters and construction


cparams = [];
cparams.eps = 1.0e-9;
pref = []; 
pref.k = 16;
narms = 0;
amp = 0.0;
% make a circle
chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref);

% make an open arc
cparams.ifclosed = 0;
chnkr2 = chunkerfunc(@(t) cos_func(t,pi,1),cparams,pref);
chnkr2 = move(chnkr2,[0;3]);

% merge arcs
chnkrs(3) = chunker();
chnkrs(1) = chnkr;
rfac = 1.1;
chnkrs(2) = move(chnkr, [0;0], [3;0], 0, rfac);
chnkrs(3) = chnkr2;
chnkrtotal = merge(chnkrs);

% form chunkgraph
cgrph = tochunkgraph(chnkrtotal);

assert(size(cgrph.verts,2) == 4)
assert(length(cgrph.echnks) == 3)
assert(cgrph.npt == chnkrtotal.npt)
assert(norm(cgrph.echnks(1).r(:) - chnkr2.r(:))<1e-14)

verts = [[0;0],[1;1],[3;0]]; edge2verts = [[1;2],[3;3]];
cgrph = chunkgraph(verts,edge2verts,{chnkr2,chnkr});
assert(isequal(cgrph.verts,verts))
assert(isequal(cgrph.edgesendverts,edge2verts))

% verify chunker shifted and scaled correctly 
vend = cgrph.echnks(1).r(:,1);
assert(norm(vend-verts(:,1)) < 1e-2);

vend = cgrph.echnks(1).r(:,end);
assert(norm(vend-verts(:,2)) < 1e-2);

end

function [r,d,d2] = cos_func(t,d,A)
% parameterization of sinusoidal boundary with period d and amplitude A
omega = 2*pi/d;
r = [t(:), A*cos(omega*t(:))].';
d = [ones(length(t),1), -omega*A*sin(omega*t(:))].';
d2 = [zeros(length(t),1), -omega^2*A*cos(omega*t(:))].';
end