function [rn,dn,d2n,dist,tn,ichn] = nearest(chnkr,ref,ich,x,u,opts)
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
%   x - precomputed Legendre nodes of order chnkr.k
%   u - the matrix created by lege.exps mapping points to coefs for order
%       chnkr.k
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

maxnewt = 200;
thresh = 1.0d-9;


if nargin < 5 || or(isempty(x),isempty(u))
    [x,~,u] = lege.exps(chnkr.k);
end
if nargin < 6
    opts = [];
end

if isfield(opts,'nitermax')
    maxnewt = opts.nitermax;
end
if isfield(opts,'thresh')
    thresh = opts.thresh;
end

dist = Inf;

if nargin < 3
    ich = 1:nch;
end

for i = 1:length(ich)
    % grab chunk data
    ii = ich(i);
    r = chnkr.rstor(:,:,ii);
    d = chnkr.dstor(:,:,ii);
    d2 = chnkr.d2stor(:,:,ii);
    h = chnkr.hstor(ii);
    rc = u*(r.');
    dc = u*(d.');
    d2c = u*(d2.');

    % starting guess on grid
    rdiff = sqrt(sum((r-ref(:)).^2,1));
    [disti,ind] = min(rdiff);
    tt = x(ind);
    ri = r(:,ind);
    di = r(:,ind);
    d2i = r(:,ind);
    
    % check the endpoints
    ts = [-1;1];
    rs = (lege.exev(ts,rc)).';
    rsdiff = sqrt(sum((rs-ref).^2,1));
    [dists,inds] = min(rsdiff);
    if dists < disti
        % continue if endpoints closer
        ri = rs(:,inds);
        if dists < dist
            dist = dists;
            ds = (lege.exev(ts,dc)).';
            d2s = (lege.exev(ts,d2c)).';
            di = ds(:,inds);
            d2i = d2s(:,inds);
            rn = ri; dn = di; d2n = d2i;
            tn = ts(inds);
            ichn = ii;
        end
        continue;
    end
        
    % otherwise run newton

    t0 = tt;
    iextra = 3;
    ifend = 0;

    for iii = 1:maxnewt
        r0 = (lege.exev(t0,rc)).';
        d0 = (lege.exev(t0,dc)).'*h;
        d20 = (lege.exev(t0,d2c)).'*h*h;

        rdiff = r0(:) - ref(:);
        dist0 = sum(rdiff.^2); 

        dprime = 2*(rdiff.'*d0(:));
        dprime2 = 2*(d0(:).'*d0(:) + d0(:).'*d20(:));
        t1 = t0 - dprime/dprime2;
        
        if t1 > 1.0
            t1 = (1.0+t0)/2;
        end
        if t1 < -1.0
            t1 = (-1.0+t0)/2;
        end
        
        err = abs(t1 - t0);

        if (err < thresh)
            ifend = ifend + 1;
        end

        t0 = t1;

        if (ifend >= iextra)
            break;
        end
    end

    if (ifend < iextra)
        warning('newton failed in nearest');
    end


    % if here, then success
    ri = (lege.exev(t1,rc)).';
    disti = sqrt(sum((ri(:)-ref(:)).^2,1));
    if disti < dist
        dist = disti;
        rn = ri;
        dn = (lege.exev(t1,dc)).'*h;
        d2n = (lege.exev(t1,d2c)).'*h*h;
        tn = t1;
        ichn = ii;
    end
end


end

