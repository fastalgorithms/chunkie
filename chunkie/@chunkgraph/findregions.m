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

    [loops] = findloops_verts(obj);

    % --- periodic geometries -------------------------------------------
    % Curves that are unbounded under the periodic identification (e.g. a
    % staircase unit cell) are handled separately from genuinely closed
    % loops. A loop's net displacement -- the sum of its edge end-minus-
    % start vectors traversed in order -- is zero for a closed loop and a
    % nonzero lattice vector (+/-dx,0)/(0,+/-dy) for a curve that only
    % closes through periodicity. The branch is gated on dx/dy being set so
    % the early findregions call during base construction (before calc_per
    % runs) falls through to the standard logic. [stage 1: single curve]
    if isa(obj_in,'chunkgraph_per') && (~isempty(obj.dx) || ~isempty(obj.dy))
        nl = numel(loops);
        isunb = false(1,nl);
        for il = 1:nl
            isunb(il) = norm(loop_displacement(obj,loops{il})) > 1e-10;
        end
        if any(isunb)
            regions = findregions_per_unbounded(obj,loops,isunb);
            return
        end
        % no unbounded loop found: fall through to the standard logic
    end

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
        %%%        imin = min(ilist);
    end

    inclusion_rels = [];

    l_out_in_old = l_out_in;

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
        %%%        imin = min(ilist);
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
%FINDREGIONS_VERTS a relatively crude method for determining the regions of 
% a chunkgraph associated with the subset of its vertices stored in iverts.
% NOTE: if iverts is not provided then all vertices will be considered. The
% Matlab routine conncomp can be used to provide subsets of vertices which 
% will define meaningful subregions.
%
% Syntax: [regions] = findloops_verts(obj, iverts);
%
% Input:
%   obj              - a chunkgraph object
%   iverts(optional) - the indices of the subset of vertices for which
%   loops
%              are to be found.
%
% Output:
%   regions - a cell array of length nloops (the number of loops
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

loops = {};
nloops = 0;

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
    
    nloops = nloops + 1;
    rcurr = {};
    rcurr = ecycle;
    loops{nloops} = rcurr;
    
    
end

end


function d = loop_displacement(obj,edges)
%LOOP_DISPLACEMENT net displacement around a loop, computed from the edge
% chunker endpoints. ~0 for a closed loop; a lattice vector for a curve
% that only closes through the periodic identification.
    d = [0;0];
    for jj = 1:numel(edges)
        e = edges(jj);
        ech = obj.echnks(abs(e));
        [r1,~] = chunkends(ech,1);
        [r2,~] = chunkends(ech,ech.nch);
        rstart = r1(:,1);
        rend   = r2(:,2);
        if e > 0
            d = d + (rend - rstart);
        else
            d = d + (rstart - rend);
        end
    end
end


function ny = mean_normal_y(obj,edges)
%MEAN_NORMAL_Y average y-component of the (orientation-adjusted) normal
% along the edges of a loop. Used to orient region 1 so its normal points
% "up" (toward the upper half-space).
    s = 0; cnt = 0;
    for jj = 1:numel(edges)
        e = edges(jj);
        nn = obj.echnks(abs(e)).n(:,:);
        nyj = nn(2,:);
        if e < 0
            nyj = -nyj;
        end
        s = s + sum(nyj);
        cnt = cnt + numel(nyj);
    end
    ny = s/cnt;
end


function regions = findregions_per_unbounded(obj,loops,isunb)
%FINDREGIONS_PER_UNBOUNDED build regions for a single open periodic curve
% that splits the plane into an upper (region 1) and lower (region 2)
% half-space. findloops_verts returns each curve in both orientations, so
% we take one representative and emit the two opposite orientations as the
% two regions (region 1 oriented so its normal points up).
    iunb = find(isunb);
    if numel(iunb) > 2
        warning(['findregions: multiple unbounded periodic curves found; ' ...
            'only a single curve is handled at this stage (layered media ' ...
            'and composites are later stages).']);
    end
    ecycle = loops{iunb(1)};
    if mean_normal_y(obj,ecycle) < 0
        ecycle = -fliplr(ecycle);     % flip so the normal points up
    end
    regions = cell(1,2);
    regions{1} = {ecycle};            % upper half-space
    regions{2} = {-fliplr(ecycle)};   % lower half-space
end