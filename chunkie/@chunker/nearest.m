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
    d = chnkr.dstor(:,:,ii);
    d2 = chnkr.d2stor(:,:,ii);
    rc = u*(r.');
    dc = u*(d.');
    d2c = u*(d2.');

    % starting guess on over sampled grid
    rover = (ainterpover*(r.')).';
    rdist0 = sqrt(sum((rover-ref(:)).^2,1));
    [disti,ind] = min(rdist0);
    tt = xover(ind);
    
    % check the endpoints
    ts = [-1;1];
    rs = (lege.exev(ts,rc)).';
    rsdist = sqrt(sum((rs-ref).^2,1));
    [dists,inds] = min(rsdist);
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
        
    % run levenberg

    t0 = tt;

    ifend = 0;

    r0 = (lege.exev(t0,rc)).';

    rdiff0 = ref(:) - r0(:);

    drc = lege.derpol(rc);
    %d2rc = lege.derpol(drc);
    thresh = thresh0;
    
    for iii = 1:maxnewt
        %d0 = (lege.exev(t0,dc)).'*h;
        d0 = (lege.exev(t0,drc)).'; % actual derivative of chunk
%        d20 = (lege.exev(t0,d2rc)).'; % actual 2nd derivative of chunk

        % newton info
%         
%         dprime = rdiff0.'*d0(:);
%         dprime2 = (d0(:).'*d0(:) + d0(:).'*d20(:));

        % levenberg instead
        
        deltat = d0(:)\(rdiff0(:));
        t1 = t0+deltat;
        
        if (t1 > 1.0)
            t1 = (1.0+t0)/2;
        end
        if (t1 < -1.0)
            t1 = (-1.0+t0)/2;
        end

        r1 = (lege.exev(t1,rc)).';
        rdiff1 = ref(:)-r1(:);

        err = min(abs(deltat),abs(sum(d0(:).*rdiff0(:))));
        if err < thresh
            ifend = ifend+1;
        end
        
        %sum(rdiff0(:).*d0(:))
        %abs(deltat)
        
        t0 = t1;
        rdiff0 = rdiff1;

%         if (err < thresh)
%             ifend = ifend + 1;
%         end
% 
        %t0 = t1;

        if (ifend >= iextra)
            break;
        end
    end
    
    % done levenberg
    
    ri = (lege.exev(t0,rc)).';
    disti = sqrt(sum((ri(:)-ref(:)).^2,1));
    if (ifend < iextra)
        warning('newton failed in nearest');
%          deltat
%          disti
%          jac
%          t0
%          rhs
%          ifend
%          sum(jac.*rhs)
    else
%        fprintf('success\n')
    end
    
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

