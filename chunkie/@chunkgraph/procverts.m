function vstruc = procverts(obj)

verts = obj.verts;
vstruc = {};

for i=1:size(verts,2)
    [inds,isgn] = vertextract(i,obj);
    vstruc{i} = {inds,isgn};
end    

end

