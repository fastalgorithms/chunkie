function cgrph = merge(cgrphs)
% merge an array of chunkgraph objects into a single chunkgraph
%
% after accumulating all verts and edgesendverts, any pair of vertices
% within 1e-14 of each other are identified, one is deleted, and
% edgesendverts is updated accordingly

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

% merge coincident vertices (within tol relative to largest vertex norm)
nverts = size(verts, 2);
tol = 1e-14 * max(vecnorm(verts));
% map(i) = canonical index for vertex i
vmap = 1:nverts;
for i = 1:nverts
    if vmap(i) ~= i
        continue   % already merged into an earlier vertex
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

% build compact vertex list and update edgesendverts
canonical = unique(vmap, 'stable');
nv_new = length(canonical);
old2new = zeros(1, nverts);
for k = 1:nv_new
    old2new(canonical(k)) = k;
end
% any vertex that was merged: redirect through its canonical representative
for i = 1:nverts
    old2new(i) = old2new(vmap(i));
end

verts = verts(:, canonical);
edgesendverts = old2new(edgesendverts);

cgrph = chunkgraph(verts, edgesendverts, fchnks);

end
