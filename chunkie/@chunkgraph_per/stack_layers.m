function cg = stack_layers(cg,merge_idx)
% merge a few chunkgraph_pers into a single chunkgraph_per
% assumes no vertices in common and no edges cross

nverts = 0;
verts = [];
edgesendverts = [];
fchnks = cell(0);
verts_per = []; 
if nargin > 1
    merge_idx = repmat(merge_idx,1,numel(cg)); 
end
for i = 1:length(cg)
    verts = [verts,cg(i).verts];
    if class(cg(i)) == "chunkgraph_per"
        edgesendverts = [edgesendverts, cg(i).edgesendverts_free+nverts];
        verts_per = [verts_per;cg(i).vert_per]; 
        merge_idx{i} = merge_idx{i}+nverts; 
    else
        edgesendverts = [edgesendverts, cg(i).edgesendverts+nverts];
    end
    nverts = nverts + size(cg(i).verts,2);
    for j = 1:size(cg(i).echnks,2)
        fchnks{end+1} = cg(i).echnks(j);
    end
end

if nargin == 2
    cg = chunkgraph_per(verts,edgesendverts,merge_idx,fchnks); 

else
    cg = chunkgraph(verts,edgesendverts,fchnks);
end

end