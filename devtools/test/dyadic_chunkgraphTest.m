dyadic_chunkgraphTest0()

function dyadic_chunkgraphTest0()
% test option to dyadically refine chunkgraph into specific corners

n = 5;
thetas = (1:n)/n*2*pi;
verts = [cos(thetas); sin(thetas)];
edge2vert = [1:n; circshift(1:n,-1)];

% vertices to refine
vert_ref = [1,4];

% init cparams
cparams = cell(1,n);
for i = 1:n
    cparams{i}.chsmall = [Inf,Inf];
end

edgeids = 1:n;
% find edges that start at a vertex that will be refined
edge_ref1 = edgeids(ismember(edge2vert(1,:),vert_ref));
for i = edge_ref1
    cparams{i}.chsmall(1) = 1e-8;
end

% find edges that end at a vertex that will be refined
edge_ref2 = edgeids(ismember(edge2vert(2,:),vert_ref));
for i = edge_ref2
    cparams{i}.chsmall(2) = 1e-8;
end

cgrph = chunkgraph(verts,edge2vert,[],cparams);

figure(1);clf
plot(cgrph,'.')

lens = chunklen(cgrph.echnks(edge_ref1(1)));
assert(lens(1) < cparams{edge_ref1(1)}.chsmall(1))
lens = chunklen(cgrph.echnks(edge_ref2(2)));
assert(lens(end) < cparams{edge_ref2(2)}.chsmall(2))

end