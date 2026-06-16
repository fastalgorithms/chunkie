function cg = merge(cgs)
%MERGE Merge an array of chunkgraph_per objects into one chunkgraph_per.
%
% authors: Tristan Goodwill, Jonathan Shaw
% merge an array of chunkgraph objects into a single chunkgraph
%
% after accumulating all verts and edgesendverts, any pair of vertices
% within 1e-14 (relative to largest vertex norm) are identified,
% edgesendverts is remapped to canonical vertices, and duplicate vertices
% are removed before calling the constructor.

% accumulate verts, edgesendverts, and edge chunkers
nverts = 0;
verts = [];
edgesendverts = [];
fchnks = cell(0);
merge_idx = []; 
for i = 1:length(cgs)
    verts = [verts, cgs(i).verts];
    edgesendverts = [edgesendverts, cgs(i).edgesendverts_free + nverts];

    if ~isempty(cgs(i).merge_idx)
         merge_idx = [merge_idx, ...
            cellfun(@(x) x + nverts, cgs(i).merge_idx, ...
            'UniformOutput', false)];
    end

    nverts = nverts + size(cgs(i).verts, 2);

    for j = 1:size(cgs(i).echnks, 2)
        fchnks{end+1} = cgs(i).echnks(j);
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
            del = []; 
            vmap(j) = i;
            if ~isempty(merge_idx) %delete indices from merge_idx
                fi = find(cellfun(@(x) (any(x==i) || any(x==j)), merge_idx), 1, 'first');
                mv = merge_idx{fi}; 
                mv(mv==i) = []; mv(mv==j) = []; 
                nmidx = numel(merge_idx); 
                for m = fi+1:nmidx
                    rv = ceil((merge_idx{m} == i) + (merge_idx{m} == j)/2); 
                    if max(rv) == 0
                        continue
                    end
                    merge_idx{m}(logical(rv)) = []; 
                    mv = [mv,merge_idx{m}]; 
                    if numel(merge_idx{m}) < 2
                        del = [del,m]; 
                    end
                end
                if ~isempty(del)
                    merge_idx(del) = []; 
                end
                merge_idx{fi} = mv; 
            end
        end
    end
end

% compact: keep only canonical vertices, re-index, and remap edgesendverts
kept = unique(vmap, 'stable');
reindex = zeros(1, nverts);
reindex(kept) = 1:numel(kept);
verts = verts(:, kept);
edgesendverts = reindex(vmap(edgesendverts));

%cleaning merge_idx:
if ~isempty(merge_idx)
   merge_idx = cellfun(@(x) reindex(x),merge_idx,'uniformoutput',false); 
end

cg = chunkgraph_per(verts,edgesendverts,merge_idx,fchnks); 

end
