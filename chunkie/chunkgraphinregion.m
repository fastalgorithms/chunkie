function [ids] = chunkgraphinregion(cg,ptsobj,opts)
%CHUNKGRAPHINREGION returns an array indicating the region number for
% each point specified by pts.
%
% Syntax: ids = chunkerinterior(cg,pts,opts)
%         ids = chunkerinterior(cg,{x,y},opts) % meshgrid version
%
% Input:
%   cg - chunkgraph object describing geometry. May also be a
%       chunkgraph_per (a periodic subclass of chunkgraph), in which case
%       the per-region interior tests are performed with the appropriate
%       period (see notes on periodicity below).
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
%  Note on periodicity (chunkgraph_per only):
%    For a chunkgraph_per, the boundary of a region need not be a closed
%    curve in the plane: it may close only up to a periodic shift. For each
%    region whose bounding edges touch merged (periodic) vertices, the
%    common period is read off from cg.vert_per and passed to
%    chunkerinterior as opts.periodic/opts.d, which then uses the
%    quasi-periodic interior test. Regions with no periodic vertices fall
%    back to the standard (non-periodic) test.
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
% see also CHUNKGRAPH, CHUNKGRAPH_PER, CHUNKERINTERIOR

% author: Travis Askham (askhamwhat@gmail.com)

if nargin < 3
    opts = [];
end

% Assign appropriate object to chnkr.
% NB: use isa(...) rather than class(...)=="chunkgraph" so that subclasses
% (e.g. chunkgraph_per) are accepted as well.
msg = "chunkgraphinregion: input 1 must be a chunkgraph (or a subclass)";
assert(isa(cg,"chunkgraph"),msg);

% periodic geometry? (chunkgraph_per carries per-vertex period info)
isper = isa(cg,"chunkgraph_per");

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

    % Empty region: skip. A periodic chunkgraph_per may have no unbounded
    % exterior loop, so findregions leaves cg.regions{1} = {} (ntot == 0).
    % There is nothing to test, and chnkrs(0) below would be an illegal
    % index. For ir == 1 this also means we do NOT apply the usual
    % "everything outside the outer boundary is region 1" default, which is
    % the correct behavior when no unbounded region exists.
    if ntot == 0
        continue
    end

    ieout = 0;
    chnkrs(ntot) = chunker(p,t,w);

    % common period for this region (NaN until a periodic vertex is seen)
    period_ir = nan;

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

            % accumulate / check the region's period from the vertices
            % this edge connects (only merged vertices have a finite period)
            if isper
                per_eid = cg.vert_per(cg.edgesendverts(:,abs(eid)));
                per_eid = per_eid(~isnan(per_eid));
                if ~isempty(per_eid)
                    assert(norm(per_eid - per_eid(1),inf) < 1e-10, ...
                        'chunkgraphinregion: unequal vertex periods on a single edge');
                    if isnan(period_ir)
                        period_ir = per_eid(1);
                    else
                        assert(abs(per_eid(1) - period_ir) < 1e-10, ...
                            'chunkgraphinregion: unequal periods within a single region');
                    end
                end
            end
        end
        
    end

    % build the options for this region's interior test, switching on the
    % periodic interior test when the region has a finite period
    opts_ir = opts;
    if isper && ~isnan(period_ir)
        if ~isstruct(opts_ir); opts_ir = struct(); end
        opts_ir.periodic = true;
        opts_ir.d = period_ir;
    end

    intmp = reshape(chunkerinterior(merge(chnkrs(1:ntot)),ptsobj,opts_ir),npts,1);

    intmp = intmp > 0;
    if ir == 1
        ids(~intmp) = ir;
    else
        ids(intmp) = ir;
    end
end
