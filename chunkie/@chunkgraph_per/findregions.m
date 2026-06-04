function [regions] = findregions(obj_in)
%FINDREGIONS_PER determines the regions of a singly-periodic chunkgraph_per.
%
% Syntax: [regions] = findregions_per(obj_in);
%
% Input:
%   obj_in - a chunkgraph_per object
%
% Output:
%   regions - a cell array of length nregions. Each region is specified by
%             a cell array of signed edge lists which traverse its boundary.
%
% Periodic geometries are represented in one canonical form:
%   layered background curves + zero or more closed cells/objects.
% Pure layered, pure closed tiling/object, and composite geometries are
% obtained by allowing either list to be empty.

msg = "findregions_per: input must be a chunkgraph_per";
assert(isa(obj_in,'chunkgraph_per'),msg);

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

regions = findregions_per_from_loops(obj,loops);

end


function regions = findregions_per_from_loops(obj,loops)
%FINDREGIONS_PER_FROM_LOOPS canonical periodic region construction.

[curves,cells] = cgper_collect_region_primitives(obj,loops);

regions = cgper_build_background_regions(curves,cells);

nbg = numel(regions);
for j = 1:numel(cells)
    regions{nbg+j} = {cells{j}};
end
end


function [curves,cells] = cgper_collect_region_primitives(obj,loops)
%CGPER_COLLECT_REGION_PRIMITIVES split unique periodic loops into
% unbounded curves and closed cells/objects.
%
% Unbounded curves have nonzero net period displacement and are oriented
% normal-up. Closed cells/objects have zero net displacement and are oriented
% counter-clockwise using cell_polygon. Both lists are ordered top to bottom
% so region numbering is stable.

loops_unique = cgper_unique_loops(loops);

curves = {};
cmy = [];

cells = {};
smy = [];

for k = 1:numel(loops_unique)
    e = loops_unique{k};
    dnet = loop_displacement(obj,e);

    if norm(dnet) > 1e-10
        if loop_normal_y(obj,e) < 0
            e = -fliplr(e);
        end

        curves{end+1} = e;
        cmy(end+1) = curve_mean_y(obj,e);
    else
        poly = cell_polygon(obj,e);
        x = poly(1,:);
        y = poly(2,:);
        A = 0.5*sum(x.*y([2:end 1]) - x([2:end 1]).*y);

        if A < 0
            e = -fliplr(e);
        end

        cells{end+1} = e;
        smy(end+1) = mean(poly(2,:));
    end
end

[~,oc] = sort(cmy,'descend');
curves = curves(oc);

[~,os] = sort(smy,'descend');
cells = cells(os);
end


function loops_unique = cgper_unique_loops(loops)
%CGPER_UNIQUE_LOOPS remove the duplicate reverse orientation of each loop.

loops_unique = {};
keys = {};

for k = 1:numel(loops)
    e = loops{k};
    if isempty(e)
        continue
    end

    key = sort(abs(e));

    isnew = true;
    for j = 1:numel(keys)
        if isequal(keys{j},key)
            isnew = false;
            break
        end
    end

    if isnew
        keys{end+1} = key;
        loops_unique{end+1} = e;
    end
end
end


function regions = cgper_build_background_regions(curves,cells)
%CGPER_BUILD_BACKGROUND_REGIONS build the layered background regions.
%
% If there are no unbounded curves, the background is the exterior of all
% closed cells. In that pure-closed case, include the reversed cell
% boundaries in region 1, preserving edge-side bookkeeping. If unbounded
% curves exist, background regions are determined only by the unbounded
% curves, and closed-cell interiors are appended as separate regions.

ncurve = numel(curves);

if ncurve == 0
    extcomps = cell(1,numel(cells));
    for j = 1:numel(cells)
        extcomps{j} = -fliplr(cells{j});
    end
    regions = {extcomps};
    return
end

regions = cell(1,ncurve+1);

regions{1} = {curves{1}};

for k = 2:ncurve
    regions{k} = {-fliplr(curves{k-1}), curves{k}};
end

regions{ncurve+1} = {-fliplr(curves{ncurve})};
end


function y = curve_mean_y(obj,edges)
%CURVE_MEAN_Y mean y-coordinate of the points making up an edge list.

s = 0;
cnt = 0;

for jj = 1:numel(edges)
    rr = obj.echnks(abs(edges(jj))).r(:,:);
    s = s + sum(rr(2,:));
    cnt = cnt + size(rr,2);
end

y = s/cnt;
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
