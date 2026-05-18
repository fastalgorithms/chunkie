function cgrph = merge(cgrphs)
% merge an array of chunkgraph_per objects into a single chunkgraph_per
%
% after accumulating all verts and edgesendverts, any pair of vertices
% within 1e-14 (relative to largest vertex norm) are identified and
% passed as merge_idx to the chunkgraph_per constructor

% accumulate verts, edgesendverts, and edge chunkers
nverts = 0;
verts = [];
edgesendverts = [];
fchnks = cell(0);
for i = 1:length(cgrphs)
    verts = [verts, cgrphs(i).verts];
    edgesendverts = [edgesendverts, cgrphs(i).edgesendverts + nverts];
    nverts = nverts + size(cgrphs(i).verts, 2);
    for j = 1:size(cgrphs(i).echnks, 2)
        fchnks{end+1} = cgrphs(i).echnks(j);
    end
end

% find coincident vertex pairs within relative tolerance
tol = 1e-14 * max(vecnorm(verts));
nverts = size(verts, 2);
merge_idx = {};
merged = false(1, nverts);
for i = 1:nverts
    if merged(i)
        continue
    end
    group = i;
    for j = i+1:nverts
        if merged(j)
            continue
        end
        if norm(verts(:,i) - verts(:,j)) < tol
            group = [group, j];
            merged(j) = true;
        end
    end
    if numel(group) > 1
        merged(i) = true;
        merge_idx{end+1} = group;
    end
end

cgrph = chunkgraph_per(verts, edgesendverts, merge_idx, fchnks);

end
