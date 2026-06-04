function [ids] = chunkgraphinregion(cg,ptsobj,opts)
%CHUNKGRAPHINREGION returns an array indicating the region number for
% each point specified by pts.
%
% Syntax: ids = chunkerinterior(cg,pts,opts)
%         ids = chunkerinterior(cg,{x,y},opts) % meshgrid version
%
% Input:
%   cg - chunkgraph or chunkgraph_per object describing geometry
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
%
% Output:
%   ids - integer array, pts(:,i) is in region ids(i)
%
% see also CHUNKGRAPH, CHUNKERINTERIOR

% author: Travis Askham (askhamwhat@gmail.com)

if nargin < 3
    opts = [];
end

msg = "chunkgraphinregion: input 1 must be chunkgraph or chunkgraph_per";
assert((class(cg) == "chunkgraph") || (class(cg) == "chunkgraph_per"),msg);

npts = get_npts(ptsobj);

ids = nan(npts,1);

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

    ntot = 0;
    for ic = 1:ncomp
        edgelist = cg.regions{ir}{ic};
        ntot = ntot + numel(edgelist);
    end

    if ntot == 0
        continue
    end

    ieout = 0;
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
    msg = "chunkgraphinregion: input 2 not a recognized type";
    error(msg);
end
end
