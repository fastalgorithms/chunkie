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
msg = "chunkgraphinregion: input 1 must be chunkgraph";
assert(class(cg) == "chunkgraph",msg);

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
    for ic = 1:ncomp
        edgelist = cg.regions{ir}{ic};
        nedge = numel(edgelist);
        chnkrs(nedge) = chunker(p,t,w);
        for ie = 1:nedge
            eid = edgelist(ie);
            if eid > 0
                if ir == 1
                    chnkrs(ie) = cg.echnks(eid);
                else
                    chnkrs(ie) = reverse(cg.echnks(eid));
                end
            else
                if ir == 1
                    chnkrs(ie) = reverse(cg.echnks(-eid));
                else
                    chnkrs(ie) = cg.echnks(-eid);
                end
            end
        end
        intmp = intmp + reshape(chunkerinterior(merge(chnkrs(1:nedge)),ptsobj,opts),npts,1);
    end

    intmp = intmp > 0;
    if ir == 1
        ids(~intmp) = ir;
    else
        ids(intmp) = ir;
    end
end
