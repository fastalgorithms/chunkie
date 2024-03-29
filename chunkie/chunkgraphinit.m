function [cgrph] = chunkgraphinit(verts,edge2verts,fchnks,cparams)

    warning('This method is deprecated. Use the chunkgraph constructor instead.');

cgrph = chunkgraph(verts,edge2verts,fchnks,cparams);

end
