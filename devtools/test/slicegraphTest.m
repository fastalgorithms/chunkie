slicegraphTest0();

function slicegraphTest0()
% create a geometry of concentric squares

% outer square
verts_out = [[1;1],[1;-1],[-1;-1],[-1;1]];

% inner square
verts_in = verts_out/2;

verts = [verts_out, verts_in];

% geometry will be a series of loops through subsets of vertices
id_vert_out = 1:4;
% call helper function
e2v_out = [id_vert_out;circshift(id_vert_out,1)];

id_vert_in = 5:8;
e2v_in = [id_vert_in;circshift(id_vert_in,1)];

% combine lists of edges and build chunkgraph
edge_2_verts = [e2v_out,e2v_in];
cgrph = chunkgraph(verts,edge_2_verts);

% slice chunkgraph and test points assigned correctly
ichs = [1,5:8];
cgrph_slc = slicegraph(cgrph, ichs);
assert(norm(cgrph_slc.r(:) - merge(cgrph.echnks(ichs)).r(:)) < 1e-14)

figure(1);
clf;
plot_regions(cgrph_slc)
axis equal


% build a test system matrix
dkern = -2*kernel('lap','d');
sysmat = chunkermat(cgrph, dkern);

% build system matrix on sliced problem
ichs = 5:8;
cgrph_slc = slicegraph(cgrph, ichs);

sysmat_slc = chunkermat(cgrph_slc, dkern);

% check that sliced matrix is correct
idslce = (cgrph.npt - cgrph_slc.npt) + (1:cgrph_slc.npt);

matdiff = norm(sysmat_slc - sysmat(idslce, idslce));
fprintf('error in sliced system matrix is %e\n',matdiff)
assert(matdiff < 1e-10)


% test getting ids of edges
ids = edgeids(cgrph,[3,4,2,1]);

% test values are correct
assert(all(sort(ids) == (1:sum([cgrph.echnks(1:4).npt]))));
% test sorted by input list
assert(~all(sort(ids) == ids));

ids = edgeids(cgrph,ichs);
% compare edgeids to manual version
assert(all(ids == idslce));

end
