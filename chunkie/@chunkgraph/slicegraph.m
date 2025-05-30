function cgrph = slicegraph(cgrph, ichs)
%SLICEGRAPH extract edges of cgrph corresponding to cgrph.echnks(ichs)

% get the edges to keep
edgesendverts = cgrph.edgesendverts(:,ichs);

% get the vertices to keep
iverts = 1:size(cgrph.verts,2);
iverts = iverts(ismember(iverts, edgesendverts));

verts = cgrph.verts(:, iverts);

% update vertex indices
edgesendvertsnew = edgesendverts;
for i = 1:size(verts,2)
    edgesendvertsnew(edgesendverts==iverts(i)) = i;
end

% assemble the sliced graph
cgrph.verts = verts;
cgrph.edgesendverts = edgesendvertsnew;
cgrph.v2emat        = build_v2emat(cgrph);
cgrph.echnks        = cgrph.echnks(ichs);
cgrph.vstruc = procverts(cgrph);
cgrph.wts = weights(cgrph);

cgrph.regions = findregions(cgrph);
end