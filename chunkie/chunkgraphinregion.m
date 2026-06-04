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
if isa(cg,"chunkgraph_per")
    if cgper_region_is_unbounded(cg)
        ids = chunkgraphinregion_per(cg,ptsobj,opts,npts);
        return
    elseif cgper_region_is_closed_tiling(cg)
        ids = chunkgraphinregion_per_closed(cg,ptsobj,opts,npts);
        return
    end
    % otherwise (closed/open, non-tiling periodic) fall through to the
    % standard closed-loop logic below
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
% the open-staircase / layered case from closed periodic geometries, which
% still go through the standard closed-loop logic.
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
    tf = norm(loop_displacement(cg,edgelist)) > 1e-10;
end


function curves = cgper_ordered_curves(cg)
%CGPER_ORDERED_CURVES recover the interface curves, ordered top to bottom,
% from the regions built by findregions. Region k (k=1..N) stores curve_k
% as its normal-up (mean normal y > 0) component.
    nreg = numel(cg.regions);
    ncurve = nreg - 1;
    curves = cell(1,ncurve);
    for k = 1:ncurve
        comp = cg.regions{k};
        pick = comp{1};
        if numel(comp) > 1
            for c = 1:numel(comp)
                if loop_normal_y(cg,comp{c}) > 0
                    pick = comp{c};
                    break
                end
            end
        end
        curves{k} = pick;
    end
end


function chnkr = cgper_build_chunker(cg,edgelist)
%CGPER_BUILD_CHUNKER merge the edge chunkers of an edge list into one
% chunker, reversing edges with negative index.
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


function ids = chunkgraphinregion_per(cg,ptsobj,opts,npts)
%CHUNKGRAPHINREGION_PER label points for one or more stacked open periodic
% curves. With N curves ordered top to bottom, a point's region is
%   region = 1 + (number of curves the point lies below),
% so region 1 is above everything and region N+1 below everything. [stage 2]
    if isempty(opts)
        optsint = struct();
    else
        optsint = opts;
    end

    curves = cgper_ordered_curves(cg);
    ncurve = numel(curves);

    countbelow = zeros(npts,1);
    for kk = 1:ncurve
        edgelist = curves{kk};
        chnkr = cgper_build_chunker(cg,edgelist);

        dnet = loop_displacement(cg,edgelist);
        optsk = optsint;
        optsk.periodic = true;
        if abs(dnet(1)) >= abs(dnet(2))
            optsk.d = cg.dx;     % x-periodic -> unbounded in y (above/below)
            vertical = true;
        else
            optsk.d = cg.dy;     % y-periodic -> unbounded in x (deferred)
            vertical = false;
        end

        inb = reshape(chunkerinterior(chnkr,ptsobj,optsk),npts,1) > 0;

        % orient: probe a point clearly below this curve so "below" maps to
        % the right logical value regardless of edge orientation / lq sign.
        rr = chnkr.r(:,:);
        if vertical
            span = max(rr(2,:)) - min(rr(2,:)) + optsk.d;
            probe = []; probe.r = [mean(rr(1,:)); min(rr(2,:)) - 10*span];
        else
            span = max(rr(1,:)) - min(rr(1,:)) + optsk.d;
            probe = []; probe.r = [min(rr(1,:)) - 10*span; mean(rr(2,:))];
        end
        inprobe = chunkerinterior(chnkr,probe,optsk) > 0;

        if inprobe
            below = inb;
        else
            below = ~inb;
        end
        countbelow = countbelow + double(below);
    end

    ids = 1 + countbelow;
end


function tf = cgper_region_is_closed_tiling(cg)
%CGPER_REGION_IS_CLOSED_TILING true if some region boundary closes only
% through periodic identification (period jumps present) with zero net
% displacement -- the case (b) "diamonds".
    tf = false;
    for ir = 1:numel(cg.regions)
        comp = cg.regions{ir};
        if ~iscell(comp); continue; end
        for c = 1:numel(comp)
            e = comp{c};
            if ~isempty(e) && norm(loop_displacement(cg,e)) < 1e-10 ...
                    && loop_max_jump(cg,e) > 1e-6
                tf = true; return
            end
        end
    end
end


function ids = chunkgraphinregion_per_closed(cg,ptsobj,opts,npts)
%CHUNKGRAPHINREGION_PER_CLOSED label points for periodic cells closed under
% tiling: region 1 is outside all cells, region k+1 is inside the k-th
% cell. Cells are disjoint, so a point is interior to at most one. [stage b]
    ids = ones(npts,1);
    if isempty(opts)
        optsint = struct();
    else
        optsint = opts;
    end
    ncell = numel(cg.regions) - 1;
    for k = 1:ncell
        edgelist = cg.regions{k+1}{1};
        chnkr = cgper_build_chunker(cg,edgelist);
        optsk = optsint;
        optsk.periodic = true;
        if ~isempty(cg.dx)
            optsk.d = cg.dx;
        else
            optsk.d = cg.dy;
        end
        inb = reshape(chunkerinterior(chnkr,ptsobj,optsk),npts,1) > 0;
        ids(inb) = k+1;
    end
end
