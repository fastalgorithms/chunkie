function [ids] = chunkgraph_perinregion(cg,ptsobj,opts)
%CHUNKGRAPH_PERINREGION returns region labels for a chunkgraph_per object.
%
% Syntax: ids = chunkgraphinregion_per(cg,pts,opts)
%         ids = chunkgraphinregion_per(cg,{x,y},opts) % meshgrid version
%
% Input:
%   cg - chunkgraph_per object describing periodic geometry
%   ptsobj - object describing the target points, can be specified as
%       * (cg.echnks(1).dim,:) array of points to test
%       * {x,y} - length 2 cell array. the points checked then have the
%           coordinates of a mesh grid [xx,yy] = meshgrid(x,y)
%       * chunker object, in which case it uses ptsobj.r(:,:)
%       * chunkgraph or chunkgraph_per object, in which case it uses
%           ptsobj.r(:,:)
%
% Optional input:
%   opts - options structure passed through to chunkerinterior
%
% Output:
%   ids - integer array, pts(:,i) is in region ids(i)
%
% Periodic region model:
%   background/layer regions = 1, ..., numel(curves)+1
%   closed-cell interiors    = numel(curves)+2, ...
%
% see also CHUNKGRAPHINREGION, CHUNKERINTERIOR
% authors: Jeremy Hoskins, Jonathan Shaw

if nargin < 3
    opts = [];
end

msg = "chunkgraphinregion_per: input 1 must be chunkgraph_per";
assert(class(cg) == "chunkgraph_per",msg);

npts = get_npts(ptsobj);

if isempty(opts)
    optsint = struct();
else
    optsint = opts;
end

[curves, closed] = cgper_collect_loops(cg);
nlayer = numel(curves);

% Layered background index.
countbelow = zeros(npts,1);
for kk = 1:nlayer
    below = cgper_below_indicator(cg,curves{kk},ptsobj,optsint,npts);
    countbelow = countbelow + double(below);
end

ids = 1 + countbelow;

% Closed sub-object interiors override the background/layer labels.
for jj = 1:numel(closed)
    inside = cgper_inside_indicator(cg,closed{jj},ptsobj,optsint,npts);
    ids(inside) = (nlayer+1) + jj;
end

end


function npts = get_npts(ptsobj)
%GET_NPTS number of target points represented by ptsobj.
if isa(ptsobj, "cell")
    assert(length(ptsobj)==2, ...
        'second input should be either 2xnpts array or length 2 cell array');
    x = ptsobj{1};
    y = ptsobj{2};
    npts = length(x)*length(y);
elseif isa(ptsobj, "chunker") || isa(ptsobj, "chunkgraph") || ...
        isa(ptsobj, "chunkgraph_per") || ...
        (isstruct(ptsobj) && isfield(ptsobj,"r"))
    npts = size(ptsobj.r(:,:),2);
elseif isnumeric(ptsobj)
    npts = size(ptsobj(:,:),2);
else
    msg = "chunkgraphinregion_per: input 2 not a recognized type";
    error(msg);
end
end

function [curves, closed] = cgper_collect_loops(cg)
%CGPER_COLLECT_LOOPS gather the distinct boundary loops from cg.regions and
% split them into unbounded periodic curves and closed cells/objects.

loops = {};
keys = {};

for ir = 1:numel(cg.regions)
    comp = cg.regions{ir};
    if ~iscell(comp)
        continue
    end

    for ic = 1:numel(comp)
        e = comp{ic};
        if isempty(e)
            continue
        end

        key = sort(abs(e));
        isnew = true;
        for g = 1:numel(keys)
            if isequal(keys{g},key)
                isnew = false;
                break
            end
        end

        if isnew
            keys{end+1} = key;
            loops{end+1} = e;
        end
    end
end

curves = {};
cmy = [];
closed = {};
smy = [];

for i = 1:numel(loops)
    e = loops{i};

    if norm(loop_displacement(cg,e)) > 1e-10
        if loop_normal_y(cg,e) < 0
            e = -fliplr(e);
        end
        curves{end+1} = e;
        cmy(end+1) = cgper_mean_y(cg,e);
    else
        poly = cell_polygon(cg,e);
        x = poly(1,:);
        y = poly(2,:);
        A = 0.5*sum(x.*y([2:end 1]) - x([2:end 1]).*y);

        if A < 0
            e = -fliplr(e);
        end

        closed{end+1} = e;
        smy(end+1) = cgper_mean_y(cg,e);
    end
end

[~,oc] = sort(cmy,'descend');
curves = curves(oc);

[~,os] = sort(smy,'descend');
closed = closed(os);
end

function below = cgper_below_indicator(cg,edgelist,ptsobj,opts,npts)
%CGPER_BELOW_INDICATOR logical array: point is below the normal-up curve.

chnkr = cgper_build_chunker(cg,edgelist);
dnet = loop_displacement(cg,edgelist);

optsk = opts;
optsk.periodic = true;

if abs(dnet(1)) >= abs(dnet(2))
    optsk.d = cg.dx;     % x-periodic -> unbounded in y
    vertical = true;
else
    optsk.d = cg.dy;     % y-periodic -> unbounded in x
    vertical = false;
end

inb = reshape(chunkerinterior(chnkr,ptsobj,optsk),npts,1) > 0;

% Orient the boolean side test by probing a point clearly below the curve.
rr = chnkr.r(:,:);
if vertical
    span = max(rr(2,:)) - min(rr(2,:)) + optsk.d;
    probe = [];
    probe.r = [mean(rr(1,:)); min(rr(2,:)) - 10*span];
else
    span = max(rr(1,:)) - min(rr(1,:)) + optsk.d;
    probe = [];
    probe.r = [min(rr(1,:)) - 10*span; mean(rr(2,:))];
end

inprobe = chunkerinterior(chnkr,probe,optsk) > 0;

if inprobe
    below = inb;
else
    below = ~inb;
end
end


function inside = cgper_inside_indicator(cg,edgelist,ptsobj,opts,npts)
%CGPER_INSIDE_INDICATOR logical array: point is inside a closed cell/object.

chnkr = cgper_build_chunker(cg,edgelist);

optsk = opts;
optsk.periodic = true;
if ~isempty(cg.dx)
    optsk.d = cg.dx;
else
    optsk.d = cg.dy;
end

inside = reshape(chunkerinterior(chnkr,ptsobj,optsk),npts,1) > 0;
end


function chnkr = cgper_build_chunker(cg,edgelist)
%CGPER_BUILD_CHUNKER merge signed edge chunkers into one chunker.

nedge = numel(edgelist);

k = cg.echnks(1).k;
t = cg.echnks(1).tstor;
w = cg.echnks(1).wstor;
p = struct("k",k);

chnkrs(nedge) = chunker(p,t,w);
for ie = 1:nedge
    eid = edgelist(ie);
    if eid > 0
        chnkrs(ie) = cg.echnks(eid);
    else
        chnkrs(ie) = reverse(cg.echnks(-eid));
    end
end

chnkr = merge(chnkrs(1:nedge));
end


function y = cgper_mean_y(cg,edgelist)
%CGPER_MEAN_Y mean y-coordinate of the points making up an edge list.

s = 0;
cnt = 0;

for jj = 1:numel(edgelist)
    rr = cg.echnks(abs(edgelist(jj))).r(:,:);
    s = s + sum(rr(2,:));
    cnt = cnt + size(rr,2);
end

y = s/cnt;
end
