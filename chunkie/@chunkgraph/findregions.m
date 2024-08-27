function [regions] = findregions(obj_in)
%FINDREGIONS determins the regions of a chunkgraph. This routine handles
% the situations involving nested chunkers inside the chunkgraph
%
% Syntax: [regions] = findregions(obj_in, iverts);
%
% Input:
%   obj_in           - a chunkgraph object
%   iverts(optional) - the indices of the subset of vertices for which regions 
%              are to be found.
%
% Output:
%   regions - a cell array of length nregions (the number of regions 
%             found). Each region is specified by a vector of 
%             indices of edges which traverse the boundary.

% author: Jeremy Hoskins
    
    obj = obj_in;
    [~, c] = find(isnan(obj.edgesendverts));
    c = unique(c);
    nnew = length(c);
    [~, nv] = size(obj.verts);
    nvnew = nv + nnew;
    verts_new = zeros(2,nvnew);
    verts_new(:,1:nv) = obj.verts;
   
    for i = 1:nnew
       verts_new(:,nv+i) = obj.echnks(c(i)).r(:,1);
       obj.edgesendverts(:, c(i)) = nv + i;
    end
    obj.verts = verts_new;
    obj.v2emat = build_v2emat(obj);
    obj.vstruc = procverts(obj);

    g = graph(obj.edgesendverts(1,:),obj.edgesendverts(2,:));
    ccomp = conncomp(g);
    
    chnkcomp = {};
    regions = {};
    
    for i=1:max(ccomp)
        inds = find(ccomp==i);
        chnkcomp{i} = inds;
        [region_comp] = findregions_verts(obj,inds);
        [region_comp] = findunbounded(obj,region_comp);
        regions{i} = region_comp;
    end
    
    gmat = zeros(numel(regions), numel(regions));
    
    for ii=1:numel(regions)
       rgna = regions{ii};
       ilist = [];
       for jj=1:numel(regions)
           if (ii ~=jj)
                rgnb = regions{jj};
                [inc] = regioninside(obj,rgnb,rgna);
                if (inc)
                    ilist = [ilist,jj];
                end
           end
           gmat(ii,ilist) = 1;
           gmat(ilist,ii) = 1;
       end
       imin = min(ilist);   
    end    
    
    ccomp_reg = conncomp(graph(gmat));
    [s,inds] = sort(ccomp_reg);
    regions = regions(inds);
    
    for ii = 1:numel(s)
        si = s(ii);
        for jj=1:(numel(s)-1)
            sj = s(jj);
            if (si == sj)
                rgna = regions{jj};
                rgnb = regions{jj+1};
                [inc] = regioninside(obj,rgna,rgnb);
                if (inc)
                    regions([jj,jj+1])= regions([jj+1,jj]);
                end
            end
        end    
    end

    
    rgns = regions;
    rgnso= {};
    
    for ii=1:max(s)
        inds = find(s==ii);
        rgnout = rgns{inds(1)};
        for jj=2:numel(inds)
            indj = inds(jj);
            [rgnout] = mergeregions(obj,rgnout,rgns{indj});
        end
        rgnso{ii} = rgnout;
    end    
    
    regions = rgnso;
    rgns = regions;
    rgnout = rgns{1};
    if (numel(rgns)>1)
        rgn2 = rgns{2};
        [rgnout] = mergeregions(obj,rgnout,rgn2);
        for ii=3:numel(rgns)
            rgn2 = rgns{ii};
            [rgnout] = mergeregions(obj,rgnout,rgn2);
        end
    end
    
    regions = rgnout;


end


function [regions] = findregions_verts(obj, iverts)
%FINDREGIONS_VERTS a relatively crude method for determining the regions of 
% a chunkgraph associated with the subset of its vertices stored in iverts.
% NOTE: if iverts is not provided then all vertices will be considered. The
% Matlab routine conncomp can be used to provide subsets of vertices which 
% will define meaningful subregions.
%
% Syntax: [regions] = findregions_verts(obj, iverts);
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

% author: Jeremy Hoskins

if (isfield(obj,'vstruc'))
    vstruc = obj.vstruc;
else
    vstruc = procverts(obj);
end    
nedge  = size(obj.edgesendverts,2);

% if the user has provided iverts, get the reduced edge list
if (nargin >1)
    v2etmp = obj.v2emat(:,iverts);
    [iinds,jinds] = find(v2etmp ~= 0);
    iinds = unique(iinds);
    edges = [iinds,-iinds];
else
    % each edge belongs to two regions (going in opposite directions)
    edges = [1:nedge,-(1:nedge)];
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
    
    if enum > 0
        iv0 = obj.edgesendverts(1,abs(enum));
        ivc = obj.edgesendverts(2,abs(enum));
    else
        iv0 = obj.edgesendverts(2,abs(enum));
        ivc = obj.edgesendverts(1,abs(enum));
    end

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
        if enum > 0
            ivc = obj.edgesendverts(2,abs(enum));
        else
            ivc = obj.edgesendverts(1,abs(enum));
        end
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

