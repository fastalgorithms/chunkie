%clear all

verts = [0,0;0,1;1,1;1,0;0.8,0.1;0.8,0.2;0.9,0.2]';
edge2verts = [-1, 1, 0, 0, 0, 0, 0; ...
               0,-1, 1, 0, 0, 0, 0; ...
               0, 0,-1, 1, 0, 0, 0; ...
               1, 0, 0,-1, 0, 0, 0; ...
              -1, 0, 1, 0, 0, 0, 0; ...
               0, 0, 0, -1, 1, 0, 0; ...
               0, 0, 0, 0,-1, 1, 0; ...
               0, 0, 0, 0, 0,-1, 1; ...
               0, 0, 0, 1, 0, 0,-1];
edge2verts = sparse(edge2verts);
verts = [verts,[3;3],[4;5],[4;2]];
edge2verts(10,8) =  -1;
edge2verts(10,9) =   1;
edge2verts(11,9) =  -1;
edge2verts(11,10) =  1;
edge2verts(12,10) = -1;
edge2verts(12,8) =   1;

fchnks    = {};
fchnks{4} = @(t) chnk.curves.fsine(t,0.2,5,0); 

prefs      = [];
prefs.chsmall = 1d-3;
[cgrph] = chunkgraphinit(verts,edge2verts,fchnks,prefs);

[inc,isgn] = vertextract(3,cgrph);
vstruc = procverts(cgrph);

rgns = findregions(cgrph);
cgrph = balance(cgrph);


