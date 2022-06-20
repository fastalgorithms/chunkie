function [cgrph] = chunkgraphinit(verts,edge2verts,echnks,prefs)
    cgrph            = chunkgraph(prefs);
    cgrph.verts      = verts;
    cgrph.edge2verts = edge2verts;
    cgrph.echnks     = echnks;
end

