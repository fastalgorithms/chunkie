function [ids] = chunkgraphinregion(cg,ptsobj,opts)
%CHUNKGRAPHINREGION returns an array indicating the region number for
% each point specified by pts.
%
% Syntax: ids = chunkerinterior(cg,pts,opts)
%         ids = chunkerinterior(cg,{x,y},opts) % meshgrid version
%
% Input:
%   cg - chunkgraph object describing geometry
%   ptsobj - object describing the target points, can be specified as
%       * (cg.echnks(1).dim,:) array of points to test
%       * {x,y} - length 2 cell array. the points checked then have the
%           coordinates of a mesh grid [xx,yy] = meshgrid(x,y)
%       * chunker object, in which case it uses ptsobj.r(:,:) 
%       * chunkgraph object, in which case it uses ptsobj.r(:,:) 
%
% Optional input:
%   opts - options structure with entries:
%       opts.fmm = boolean, use FMM 
%       opts.flam = boolean, use FLAM routines
%       opts.axissym = boolean, chunker is axissymmetric
%  Note on the default behavior: 
%    by default it tries to use the fmm if it exists, if it doesn't
%    then unless explicitly set to false, it tries to use flam
%
% Output:
%   ids - integer array, pts(:,i) is in region ids(i)
%
% Examples:
%   verts = [1 0 -1 0; 0 1 0 -1]; edgesendverts = [1:3, 3, 4; 2:3, 1, 4, 1];
%   cg = chunkgraph(verts,edgesendverts);
%   pts = 2*randn(2,100);
%   ids = chunkgraphinregion(cg,pts);
%
% see also CHUNKGRAPH, CHUNKERINTERIOR

% author: Travis Askham (askhamwhat@gmail.com)

if nargin < 3
    opts = [];
end

% Assign appropriate object to chnkr
msg = "chunkgraphinregion: input 1 must be chunkgraph or chunkgraph_per";
assert((class(cg) == "chunkgraph") || (class(cg) == "chunkgraph_per"),msg);

% Figure out size of ids array based on ptsobj
if isa(ptsobj, "cell")
    assert(length(ptsobj)==2,'second input should be either 2xnpts array or length 2 cell array');
    x = ptsobj{1};
    y = ptsobj{2};
    npts = length(x)*length(y);
elseif isa(ptsobj, "chunker") || isa(ptsobj, "chunkgraph") || ...
        (isstruct(ptsobj) && isfield(ptsobj,"r"))
    npts = size(ptsobj.r(:,:),2);
elseif isnumeric(ptsobj)
    npts = size(ptsobj(:,:),2);
else
    msg = "chunkgraphinregion: input 2 not a recognized type";
    error(msg);
end

ids = nan(npts,1);

% periodic, unbounded geometries (e.g. a single open staircase cell) are
% labeled by the periodic interior test rather than the closed-loop logic
% below. [stage 1: single open curve -> region 1 above, region 2 below]
if isa(cg,"chunkgraph_per") && cgper_region_is_unbounded(cg)
    ids = chunkgraphinregion_per(cg,ptsobj,opts,npts);
    return
end

% loop over regions and use chunkerinterior to label
% TODO: make a more efficient version 

nedge0 = numel(cg.echnks);
if nedge0 > 0
    k = cg.echnks(1).k;
    t = cg.echnks(1).tstor;
    w = cg.echnks(1).wstor;
    p = struct("k",k);
else
    ids = [];
    return
end

nr = numel(cg.regions);
for ir = 1:nr
    ncomp = numel(cg.regions{ir});
    intmp = zeros(npts,1);
    
    ntot = 0;
    for ic = 1:ncomp
            edgelist = cg.regions{ir}{ic};
        nedge = numel(edgelist);
        ntot = ntot + nedge;
    end
    ieout = 0;
    if ntot == 0 %PATCH... CONSIDER FIXING LATER
        continue
    end
    chnkrs(ntot) = chunker(p,t,w);
    
    for ic = 1:ncomp
        edgelist = cg.regions{ir}{ic};
        nedge = numel(edgelist);

        for ie = 1:nedge
            ieout = ieout + 1;
            eid = edgelist(ie);
            if eid > 0
                if ir == 1
                    chnkrs(ieout) = cg.echnks(eid);
                else
                    chnkrs(ieout) = reverse(cg.echnks(eid));
                end
            else
                if ir == 1
                    chnkrs(ieout) = reverse(cg.echnks(-eid));
                else
                    chnkrs(ieout) = cg.echnks(-eid);
                end
            end
        end
        
    end
    intmp = reshape(chunkerinterior(merge(chnkrs(1:ntot)),ptsobj,opts),npts,1);

    intmp = intmp > 0;
    if ir == 1
        ids(~intmp) = ir;
    else
        ids(intmp) = ir;
    end
end

end


function tf = cgper_region_is_unbounded(cg)
%CGPER_REGION_IS_UNBOUNDED true if region 1's boundary is an unbounded
% periodic curve (nonzero net displacement around the loop). Distinguishes
% the open-staircase case from closed periodic geometries, which still go
% through the standard closed-loop logic.
    tf = false;
    if isempty(cg.regions) || isempty(cg.regions{1})
        return
    end
    comp = cg.regions{1};
    if ~iscell(comp)
        return
    end
    edgelist = comp{1};
    if isempty(edgelist)
        return
    end
    tf = norm(cgper_net_disp(cg,edgelist)) > 1e-10;
end


function dnet = cgper_net_disp(cg,edgelist)
%CGPER_NET_DISP net displacement around a boundary edge list, from the edge
% chunker endpoints. (Mirrors loop_displacement in findregions; factor into
% a shared method when stages (b)/(f) land.)
    dnet = [0;0];
    for ie = 1:numel(edgelist)
        eid = edgelist(ie);
        ech = cg.echnks(abs(eid));
        [r1,~] = chunkends(ech,1);
        [r2,~] = chunkends(ech,ech.nch);
        if eid > 0
            dnet = dnet + (r2(:,2) - r1(:,1));
        else
            dnet = dnet + (r1(:,1) - r2(:,2));
        end
    end
end


function ids = chunkgraphinregion_per(cg,ptsobj,opts,npts)
%CHUNKGRAPHINREGION_PER label points for a single open periodic curve.
% Region 1 is the upper half-space, region 2 the lower. [stage 1]
    ids = nan(npts,1);

    if isempty(opts)
        optsint = struct();
    else
        optsint = opts;
    end

    % reconstruct the unit-cell curve from region 1's boundary edges
    edgelist = cg.regions{1}{1};
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

    % which axis is periodic? -> sets the period and the unbounded direction
    dnet = cgper_net_disp(cg,edgelist);
    optsint.periodic = true;
    if abs(dnet(1)) >= abs(dnet(2))
        optsint.d = cg.dx;     % x-periodic -> unbounded in y (above/below)
        vertical = true;
    else
        optsint.d = cg.dy;     % y-periodic -> unbounded in x (deferred)
        vertical = false;
    end

    % split the targets with the periodic interior test
    inb = reshape(chunkerinterior(chnkr,ptsobj,optsint),npts,1) > 0;

    % orient labels by probing a point on the region-1 (upper) side. This
    % keeps the labeling independent of edge orientation / the lq sign.
    rr = chnkr.r(:,:);
    if vertical
        span = max(rr(2,:)) - min(rr(2,:)) + optsint.d;
        probe = []; probe.r = [mean(rr(1,:)); max(rr(2,:)) + 10*span];
    else
        span = max(rr(1,:)) - min(rr(1,:)) + optsint.d;
        probe = []; probe.r = [max(rr(1,:)) + 10*span; mean(rr(2,:))];
    end
    inprobe = chunkerinterior(chnkr,probe,optsint) > 0;

    if inprobe
        ids(inb)  = 1;
        ids(~inb) = 2;
    else
        ids(~inb) = 1;
        ids(inb)  = 2;
    end
end
