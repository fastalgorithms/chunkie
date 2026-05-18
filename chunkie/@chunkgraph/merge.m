function cgrph = merge(cgrphs)
% merge an array of chunkgraph objects into a single chunkgraph
%
% after accumulating all verts and edgesendverts, any pair of vertices
% within 1e-14 (relative to largest vertex norm) are identified, and
% edgesendverts is updated so all references point to one canonical vertex
% per coincident group. the full vertex list is passed to the constructor
% so that edge snapping works correctly.

% accumulate verts, edgesendverts, and edge chunkers
nverts = 0;
verts = [];
edgesendverts = [];
fchnks = cell(0);
for i = 1:length(cgrphs)
    verts = [verts, cgrphs(i).verts];
    if isa(cgrphs(i), 'chunkgraph_per')
        edgesendverts = [edgesendverts, cgrphs(i).edgesendverts_free + nverts];
    else
        edgesendverts = [edgesendverts, cgrphs(i).edgesendverts + nverts];
    end
    nverts = nverts + size(cgrphs(i).verts, 2);
    for j = 1:size(cgrphs(i).echnks, 2)
        fchnks{end+1} = cgrphs(i).echnks(j);
    end
end

% find coincident vertex groups within relative tolerance
tol = 1e-14 * max(vecnorm(verts));
nverts = size(verts, 2);
vmap = 1:nverts;
for i = 1:nverts
    if vmap(i) ~= i
        continue
    end
    for j = i+1:nverts
        if vmap(j) ~= j
            continue
        end
        if norm(verts(:,i) - verts(:,j)) < tol
            vmap(j) = i;
        end
    end
end

% remap edgesendverts to canonical vertices; keep full vertex list
edgesendverts = vmap(edgesendverts);

cgrph = chunkgraph(verts, edgesendverts, fchnks);

end
