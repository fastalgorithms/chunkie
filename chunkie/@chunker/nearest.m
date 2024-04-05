function [rn,dn,d2n,dist,tn,ichn] = nearest(chnkr,ref,ich,opts,u, ...
    xover,ainterpover)
%NEAREST Find nearest point on chunker for given reference point. 
% If ich provided, only checks chunks in ich. If not, checks all chunks.
%
% Syntax: [rn,dn,d2n,dist,tn,ichn] = nearest(chnkr,ref,ich,x,u,opts)
%
% Input:
%   chnkr - chunker object
%   ref - the reference point whose distance to curve is being computed
%
% Optional input:
%   ich - vector of chunks to check (instead of checking all)
%   u - the matrix created by lege.exps mapping point values to coefs for 
%       order chnkr.k
%   xover - precomputed oversampled (order 2*chnkr.k) lege grid
%   ainterpover - precomputed interpolation matrix to oversampled lege grid
%   opts - options structure
%       opts.thresh - threshold for newton (1e-9)
%       opts.nitermax - maximum iterations for newton (200)
%
% Output:
%   rn - nearest point
%   dn - derivative at nearest point
%   d2n - second derivative at nearest point
%   dist - distance
%   tn - point in [-1,1] on ichn
%   ichn - chunk containing nearest point
%
% Examples:
%   % check all chunks for nearest point
%   [rn,~,~,dist] = nearest(chnkr,ref);
%   % only check some chunks
%   ich = 3:10;
%   rn = nearest(chnkr,ref,ich);
%

% author: Travis Askham (askhamwhat@gmail.com)

maxnewt = 15;
thresh0 = 1.0d-9;
iextra = 3;
nch = chnkr.nch;

if nargin < 5 || isempty(u)
    [~,~,u] = lege.exps(chnkr.k);
end
if nargin < 7 || or(isempty(xover),isempty(ainterpover))
    xover = lege.exps(2*chnkr.k);
    ainterpover = lege.matrin(chnkr.k,xover);
end
if nargin < 4
    opts = [];
end

if isfield(opts,'nitermax')
    maxnewt = opts.nitermax;
end
if isfield(opts,'thresh')
    thresh0 = opts.thresh;
end

dist = Inf;

if nargin < 3
    ich = 1:nch;
end

for i = 1:length(ich)
    % grab chunk data
    ii = ich(i);
    r = chnkr.rstor(:,:,ii);

    [ti,ri,di,d2i,disti] = chnk.chunk_nearparam(r,ref,opts,chnkr.tstor,u);
    disti = sqrt(disti);
    if disti < dist
        dist = disti;
        rn = ri;
        dn = di;
        d2n = d2i;
        tn = ti;
        ichn = ii;
    end
end


end

