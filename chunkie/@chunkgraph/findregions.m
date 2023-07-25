function [regions] = findregions(obj,iverts)
%FINDREGIONS a relatively crude method for determining the regions of 
% a chunkgraph associated with the subset of its vertices stored in iverts.
% NOTE: if iverts is not provided then all vertices will be considered. The
% Matlab routine conncomp can be used to provide subsets of vertices which 
% will define meaningful subregions.
%
% Syntax: [regions] = findregions(obj,iverts);
%
% Input:
%   obj              - a chunkgraph object
%   iverts(optional) - the indices of the subset of vertices for which regions 
%              are to be found.
%
% Output:
%   regions - a cell array of length nregions (the number of regions 
%             found). Each region is specified by a vector of 
%             indices of edges which traverse the boundary.
%  
%
%

% author: Jeremy Hoskins

if (isfield(obj,'vstruc'))
    vstruc = obj.vstruc;
else
    vstruc = procverts(obj);
end    
nedge  = size(obj.verts2edge,1);
e2v    = obj.verts2edge;

% each edge belongs to two regions (going in opposite directions)
edges = [1:nedge,-(1:nedge)];

% if the user has provided iverts, get the reduced edge list
if (nargin >1)
    e2vtmp = e2v(:,iverts);
    [iinds,jinds] = find(e2vtmp ~= 0);
    iinds = unique(iinds);
    edges = [iinds,-iinds];
end    

regions = {};
nregions = 0;

% Regions are obtained by picking an edge (including orientation) and 
% constructing a path by choosing the next edge (counterclockwise) at 
% each subsequent vertex. This will give a region with no edges passing 
% through it (unless the graph isn't planar...). The edges are then deleted
% from the stack.
while (numel(edges)>0)
    
    enum   = edges(1);
    edges(1) = [];
    estart = enum;
    ecycle = [enum];
    
    iv0 = find(e2v(abs(enum),:)==-sign(enum));
    ivc = find(e2v(abs(enum),:)== sign(enum));
    ifdone = false;

    while (~ifdone)
        inds = find(vstruc{ivc}{1} == abs(enum));
        irel = find(vstruc{ivc}{2}(inds) == sign(enum));
        ind  = inds(irel);
        if (ind < numel(vstruc{ivc}{1}))
            enext = vstruc{ivc}{1}(ind+1);
            esign = vstruc{ivc}{2}(ind+1);
        else
            enext = vstruc{ivc}{1}(1);
            esign = vstruc{ivc}{2}(1);
        end
        enum = -esign*enext;
        ivc = find(e2v(abs(enum),:)== sign(enum));
        if (enum == estart)
            ifdone = true;
        else
            ecycle = [ecycle,enum];
            edges(edges==enum) = [];
        end  
    end  
    
    nregions = nregions + 1;
    rcurr = {};
    rcurr{1} = ecycle;
    regions{nregions} = rcurr;
    
    
end

end

