function pseudo_normals = get_pseudo_normals(verts)
    [~, nv] = size(verts);
    verts_ext = [verts, verts(:,1)];
    d = verts_ext(:,2:end) - verts_ext(:,1:nv);
    
end