function edge_regs = find_edge_regions(obj)
%FIND_EDGE_REGIONS find regions on either side of the edges
% in a chunkgraph.
%
% Syntax: edge_regs = find_edge_regions(obj)
%
% Input:
%   obj - chunkgraph object describing curve
%   
% Output:
%   edge_regs - (2,nedge) array where nedges
%     is the number of edges in the chunkgraph. 
%     edge_reg(1,i) is index of the region in which the 
%     normal to edge i is pointing and edge_reg(2,i)
%     is the index corresponding to the negative of
%     the normal.
%
% author: Tristan Goodwill


    edge_regs = zeros(2,size(obj.edgesendverts,2));
    nreg = length(obj.regions);
    for i = 1:nreg
        id_edge = obj.regions{i}{1};
        edge_regs(1,id_edge(id_edge>0)) = i;
        edge_regs(2,-id_edge(id_edge<0)) = i;
    end
end
