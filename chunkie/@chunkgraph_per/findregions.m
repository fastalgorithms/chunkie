function [regions] = findregions(obj_in)
%FINDREGIONS in @chunkgraph_per determines the regions of chunkgraph_per object.
%
% Syntax: [regions] = findregions(obj_in);
%
% Input:
%   obj_in - a chunkgraph_per object
%
% Output:
%   regions - a cell array of length nregions. Each region is specified by
%             a cell array of signed edge lists which traverse its boundary.
%
% Periodic geometries are represented as
%   layered background curves + closed objects.

%authors: Jeremy Hoskins, Jonathan Shaw 

msg = "findregions_per: input must be a chunkgraph_per";
assert(isa(obj_in,'chunkgraph_per'),msg);

%regions for open components: 
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
[unbnd,bnd] = classify_loops(obj,loops);
regions = make_bg_regions(unbnd,bnd); 

nbg = numel(regions);
for j = 1:numel(bnd)
    regions{nbg+j} = {bnd{j}};
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
function [unbnd,bnd] = classify_loops(obj,loops)
%CLASSIFY_LOOPS sort unique loops into
% unbounded curves and bounded+closed objects.

loops_unique = get_unique_loops(loops);

unbnd = {};
cmy = [];

bnd = {};
smy = [];
dtol = 1e-10; 

for k = 1:numel(loops_unique)
    e = loops_unique{k};
    dnet = loop_displacement(obj,e);

    if norm(dnet) > dtol
        if loop_normal_y(obj,e) < 0
            e = -fliplr(e);
        end

        unbnd{end+1} = e;
        cmy(end+1) = loop_mean_y(obj,e);
    else
        poly = cell_polygon(obj,e);
        x = poly(1,:);
        y = poly(2,:);
        A = 0.5*sum(x.*y([2:end 1]) - x([2:end 1]).*y); %signed area

        if A < 0
            e = -fliplr(e);
        end

        bnd{end+1} = e;
        smy(end+1) = mean(y);
    end
end

%sort by mean y-component:
[~,oc] = sort(cmy,'descend');
unbnd = unbnd(oc);

[~,os] = sort(smy,'descend');
bnd = bnd(os);
end

function loops_unique = get_unique_loops(loops)
%GET_UNIQUE_LOOPS remove the duplicate reverse orientation of each loop.

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
%}
function regions = make_bg_regions(unbnd,bnd)
%MAKE_BG_REGIONS build the layered background regions.
%
% If there are no unbounded curves, the background is the exterior of all
% closed objects. 

nunbnd = numel(unbnd);

if nunbnd == 0
    extcomps = cell(1,numel(bnd));
    for j = 1:numel(bnd)
        extcomps{j} = -fliplr(bnd{j});
    end
    regions = {extcomps};
    return
end

regions = cell(1,nunbnd+1);

regions{1} = {unbnd{1}};

for k = 2:nunbnd
    regions{k} = {-fliplr(unbnd{k-1}), unbnd{k}};
end

regions{nunbnd+1} = {-fliplr(unbnd{nunbnd})};
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
