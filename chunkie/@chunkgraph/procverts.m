function vstruc = procverts(obj)
%PROCVERTS for each vertex in a chunkgraph object, sorts the edges
% connecting to that vertex in counterclockwise order. The results are 
% stored in a struct.
%
% Syntax: vstruc = procverts(cgrph);
%
% Input:
%   cgrph    - chunkgraph object
%
% Output:
%   vstruc   - a cell array with one cell per each vertex in the 
%              chunkgrph. Each cell consists of a list of edge numbers 
%              which terminate at that vertex, and a vector of signs 
%              which are positive if the corresponding edge ends at the 
%              vertex and negative if the edge begins at the vertex.
%  
%
%

% author: Jeremy Hoskins

verts = obj.verts;
vstruc = {};

% loop over all the vertices and call vertextract
for i=1:size(verts,2)
    [inds,isgn] = vertextract(i,obj);
    vstruc{i} = {inds,isgn};
end    

end

