function [ids] = chunkgraphinregion(chnkobj,ptsobj,opts)
%CHUNKERINTERIOR returns an array indicating whether each point specified 
% % by pts is inside the domain. Assumes the domain is closed.
%
% Syntax: in = chunkerinterior(chnkobj,pts,opts)
%         in = chunkerinterior(chnkobj,{x,y},opts) % meshgrid version
%
% Input:
%   chnkobj - chunker object or chunkgraph object describing geometry
%   ptsobj - object describing the target points, can be specified as
%       * (chnkr.dim,:) array of points to test
%       * {x,y} - length 2 cell array. the points checked then have the
%           coordinates of a mesh grid [xx,yy] = meshgrid(x,y)
%       * chunker object, in which case it uses chunker.r(:,:) 
%       * chunkgraph object, in which case it uses chunkgraph.r(:,:) 
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
%   in - logical array, if in(i) is true, then pts(:,i) is inside the
%       domain or for a mesh grid [xx(i); yy(i)] is inside the domain.
%
% Examples:
%   chnkr = chunkerfunc(@(t) starfish(t));
%   pts = 2*randn(2,100);
%   in = chunkerinterior(chnkr,pts);
%

% author: Travis Askham (askhamwhat@gmail.com)

if nargin < 3
    opts = [];
end

% Assign appropriate object to chnkr
msg = "chunkgraphinregion: input 1 must be chunkgraph";
assert(class(chnkobj) == "chunkgraph",msg);

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

nedge0 = numel(chnkobj.echnks);
if nedge0 > 0
    k = chnkobj.echnks(1).k;
    t = chnkobj.echnks(1).tstor;
    w = chnkobj.echnks(1).wstor;
    p = struct("k",k);
else
    ids = [];
    return
end

nr = numel(chnkobj.regions);
for ir = 1:nr
    ncomp = numel(chnkobj.regions{ir});
    intmp = zeros(npts,1);
    for ic = 1:ncomp
        edgelist = chnkobj.regions{ir}{ic};
        nedge = numel(edgelist);
        chnkrs(nedge) = chunker(p,t,w);
        for ie = 1:nedge
            eid = edgelist(ie);
            if eid > 0
                if ir == 1
                    chnkrs(ie) = chnkobj.echnks(eid);
                else
                    chnkrs(ie) = reverse(chnkobj.echnks(eid));
                end
            else
                if ir == 1
                    chnkrs(ie) = reverse(chnkobj.echnks(-eid));
                else
                    chnkrs(ie) = chnkobj.echnks(-eid);
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
