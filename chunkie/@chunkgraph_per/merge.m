function cgrph = merge(cgrphs)
%MERGE Merge an array of chunkgraph_per objects into one chunkgraph_per.
%
% This preserves existing periodic vertex identifications from each input
% chunkgraph_per and also identifies coincident vertices between different
% input graphs.

% accumulate verts, free edge endpoints, edge chunkers, and periodic merges
nverts = 0;
verts = [];
edgesendverts = [];
fchnks = cell(0);
merge_idx = {};

for i = 1:length(cgrphs)

    % Accumulate vertices.
    verts = [verts, cgrphs(i).verts];

    % For chunkgraph_per inputs, use the unmerged/free edge endpoints.
    % The constructor will re-apply all periodic identifications.
    if isa(cgrphs(i), 'chunkgraph_per') && ~isempty(cgrphs(i).edgesendverts_free)
        edgesendverts = [edgesendverts, cgrphs(i).edgesendverts_free + nverts];
    else
        edgesendverts = [edgesendverts, cgrphs(i).edgesendverts + nverts];
    end

    % Preserve periodic merge groups already present in each object.
    if isa(cgrphs(i), 'chunkgraph_per') && ...
            ~isempty(cgrphs(i).merge_idx) && ...
            ~isempty(cgrphs(i).merge_idx{1})
        merge_idx = [merge_idx, ...
            cellfun(@(x) x + nverts, cgrphs(i).merge_idx, ...
            'UniformOutput', false)];
    end

    % Accumulate edge chunkers.
    for j = 1:size(cgrphs(i).echnks, 2)
        fchnks{end+1} = cgrphs(i).echnks(j);
    end

    nverts = nverts + size(cgrphs(i).verts, 2);
end

% Find additional coincident vertices between/among the accumulated graphs.
% These are ordinary geometric coincidences, not periodic translations.
if ~isempty(verts)
    tol = 1e-14 * max(1, max(vecnorm(verts)));
else
    tol = 0;
end

nverts = size(verts, 2);
used = false(1, nverts);

for i = 1:nverts
    if used(i)
        continue
    end

    group = i;
    for j = i+1:nverts
        if used(j)
            continue
        end

        if norm(verts(:,i) - verts(:,j)) < tol
            group = [group, j];
            used(j) = true;
        end
    end

    if numel(group) > 1
        used(group) = true;
        merge_idx{end+1} = group;
    end
end

% Constructor expects a nonempty cell array whose first cell can be empty
% when there are no periodic/coincident merges.
if isempty(merge_idx)
    merge_idx = {[]};
end

cgrph = chunkgraph_per(verts, edgesendverts, merge_idx, fchnks);

end
