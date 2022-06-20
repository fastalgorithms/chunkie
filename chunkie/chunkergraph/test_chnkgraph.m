clear all

verts = [0,0;0,1;1,1;1,0]';
edge2verts = [-1,1,0,0;0,-1,1,0;0,0,-1,1;1,0,0,-1;-1,0,1,0];
echnks     = chunker.empty;
prefs      = [];
[cgrph] = chunkgraphinit(verts,edge2verts,echnks,prefs);
