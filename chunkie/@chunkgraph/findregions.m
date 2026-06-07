function [regions] = findregions(obj_in)
%FINDREGIONS determines the regions of a chunkgraph_per object. 
%
% Syntax: [regions] = findregions(obj_in, iverts);
%
% Input:
%   obj_in           - a chunkgraph
%   iverts(optional) - the indices of the subset of vertices for which
%              regions are to be found.
%
% Output:
%   regions - a cell array of length nregions. Each region is specified by
%             a cell array of signed edge lists which traverse its boundary.

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
[loops] = findloops_verts(obj);
[li,lo] = findunbounded_loop(obj,loops);

%%%%
%%%%        .   .   .   check which lo loops are in li
%%%%

l_out_in = cell(1,numel(li));
l_out_ifin = ones(numel(lo),1);

for ii=1:numel(li)
    li_a = li{ii};
    ilist = [];
    for jj=1:numel(lo)
        lo_b = lo{jj};
        [inc] = loopinside(obj,li_a,lo_b);
        if (inc)
            ilist = [ilist,jj];
            l_out_ifin(jj) = 0;
        end
    end
    l_out_in{ii} = ilist;
end

inclusion_rels = [];

for ii=1:numel(lo)
    lo_a = lo{ii};
    for jj=1:numel(lo)
        if (ii ~=jj)
            lo_b = lo{jj};
            [inc] = loopinside(obj,lo_a,lo_b);
            if (inc)
                inclusion_rels = [inclusion_rels, [ii;jj]];
            end
        end
    end
end

for ii=1:numel(l_out_in)
    ilist = l_out_in{ii};
    dels = zeros(size(ilist));
    for kk=1:size(inclusion_rels,2)
        i_up = inclusion_rels(1,kk);
        i_dw = inclusion_rels(2,kk);
        iind_up = find(ilist == i_up);
        if (numel(iind_up) ~= 0)
            iind_dw = find(ilist == i_dw);
            if (numel(iind_dw)~=0)
                dels(iind_dw) = dels(iind_dw) + 1;
            end
        end
    end
    inds = find(dels);
    ilist(inds) = [];
    l_out_in{ii} = ilist;
end

regions = {};
for ii = 1:numel(li)
    rg = [li(ii),lo(l_out_in{ii})];
    regions{ii+1} = rg(:).';
end

inds_unbound = find(l_out_ifin);
regions{1} = lo(inds_unbound);

end

function [loops] = findloops_verts(obj, iverts)
%FINDLOOPS_VERTS a relatively crude method for determining the loops of a
% chunkgraph associated with the subset of its vertices stored in iverts.
%
% NOTE: if iverts is not provided then all vertices are considered. The
% Matlab routine conncomp can be used to provide subsets of vertices which
% define meaningful subregions.

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

loops = {};
nloops = 0;

% Regions are obtained by picking an edge (including orientation) and
% constructing a path by choosing the next edge (counterclockwise) at each
% subsequent vertex. The edges are then deleted from the stack.
while (numel(edges)>0)

    enum   = edges(1);
    edges(1) = [];
    estart = enum;
    ecycle = [enum];

    if enum > 0
        ivc = obj.edgesendverts(2,abs(enum));
    else
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

    nloops = nloops + 1;
    loops{nloops} = ecycle;
end
end
