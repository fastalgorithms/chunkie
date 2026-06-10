function [ids] = chunkgraph_perinregion(obj,ptsobj,opts)
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
%   closed-object interiors    = numel(curves)+2, ...
%
% see also CHUNKGRAPHINREGION, CHUNKERINTERIOR
% authors: Jeremy Hoskins, Jonathan Shaw

if nargin < 3
    opts = []; 
end

msg = "chunkgraphinregion_per: input 1 must be chunkgraph_per";
assert(class(obj) == "chunkgraph_per",msg);

npts = get_npts(ptsobj);

[unbnd, bnd] = classify_reg_loops(obj); %id loops as unbounded curves or closed+bounded objects
nlayer = numel(unbnd);
opts.periodic = true; 
opts.dx = obj.dx;

%region syntax: increase region # in downward direction for layered media:
ids = ones(npts,1);
for l = 1:nlayer
    chnkr = make_chunker(obj,unbnd{l});
    iidx = reshape(chunkerinterior(chnkr,ptsobj,opts),npts,1) > 0;
    ids = ids + iidx;
end

%closed object interiors override the background/layer labels:
for b = 1:numel(bnd)
    chnkr = make_chunker(obj,bnd{b}); 
    iidx = reshape(chunkerinterior(chnkr,ptsobj,opts),npts,1) > 0;
    ids(iidx) = (nlayer+1) + b;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

function [unbnd, bnd] = classify_reg_loops(cg)
%CLASSIFY_REG_LOOPS gather the distinct boundary loops from cg.regions and
% split them into unbounded periodic curves and bounded+closed objects.

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

unbnd = {};
cmy = [];
bnd = {};
smy = [];

for i = 1:numel(loops)
    e = loops{i};

    if norm(loop_displacement(cg,e)) > 1e-10
        if loop_normal_y(cg,e) < 0
            e = -fliplr(e);
        end
        unbnd{end+1} = e;
        cmy(end+1) = loop_mean_y(cg,e);
    else
        poly = cell_polygon(cg,e);
        x = poly(1,:);
        y = poly(2,:);
        A = 0.5*sum(x.*y([2:end 1]) - x([2:end 1]).*y);

        if A < 0
            e = -fliplr(e);
        end

        bnd{end+1} = e;
        smy(end+1) = loop_mean_y(cg,e);
    end
end

[~,oc] = sort(cmy,'descend');
unbnd = unbnd(oc);

[~,os] = sort(smy,'descend');
bnd = bnd(os);
end

function chnkr = make_chunker(cg,edgelist)
%MAKE_CHUNKER merge signed edge chunkers into one chunker.

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

