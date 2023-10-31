function [ts,rs,ds,d2s,dist2s] = chunk_nearparam(rval,pts,opts,t,u)
%CHNK.CHUNK_NEARPARAM Find nearest point on a single chunk for given 
% reference point. This routine is not intended to be user-callable
%
% Syntax: [rn,dn,d2n,dist,tn,ichn] = nearest(rval,pts,opts,t,u)
%
% Input:
%   rval - dim x k array of chunk positions
%   pts - the reference point whose distance to curve is being computed
%
% Optional input:
%   u - the matrix created by lege.exps mapping point values to coefs for 
%       order chnkr.k
%   opts - options structure
%       opts.thresh - threshold for newton (1e-9)
%       opts.nitermax - maximum iterations for newton (200)
%
% Output:
%   ts - ts(j) in [-1,1] is preimage on chunk that is closest to pts(:,j)
%   rs - rs(:,j) chunk coordinates at ts(j)
%

% author: Travis Askham (askhamwhat@gmail.com)

maxnewt = 15;
thresh0 = 1.0d-14;

if nargin < 4 || isempty(t) || nargin < 5 || isempty(u)
    [t,~,u] = lege.exps(chnkr.k);
end

if nargin < 3
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

dim = size(pts,1);
npts = numel(pts)/dim;
pts = reshape(pts,dim,npts);

k = size(rval,2);

rc = u*(rval.');
drc = lege.derpol(rc); drc = [drc;zeros(1,dim)];
d2rc = lege.derpol(drc); d2rc = [d2rc;zeros(1,dim)];
cfs = [rc.';drc.';d2rc.'].';

dist2all = reshape(sum( abs(reshape(pts,dim,1,npts) ...
    - reshape(rval,dim,k,1)).^2, 1),k,npts);
[dist2all,ipt] = min(dist2all,[],1);

rs = zeros(dim,npts);
ds = zeros(dim,npts);
d2s = zeros(dim,npts);
dist2s = zeros(npts,1);
ts = zeros(npts,1);

thresh = thresh0 * (k^2*sum(abs(drc(:))) + k*sum(abs(rc(:))));
   
for i = 1:npts
    
    ref = pts(:,i);

    % closest point on grid

    t0 = t(ipt(i));

    all0 = lege.exev(t0,cfs);
    r0 = all0(1:dim);
    d0 = all0(dim+1:2*dim);
    d20 = all0(2*dim+1:end);

    ts(i) = t(ipt(i));
    rs(:,i) = r0;
    ds(:,i) = d0;
    d2s(:,i) = d20;
    dist2s(i) = dist2all(i);

    % try newton 

    rdiff0 = r0(:) - ref;
    dprime = rdiff0.'*d0(:);
    dprime2 = d0(:).'*d0(:) + d0(:).'*d20(:);

    newtsuccess = false;

    for iii = 1:maxnewt

        % Newton step
        
        deltat = - dprime/dprime2;

        t0 = t0+deltat;
        
        if (t0 > 1.0)
            t0 = 1;
        end
        if (t0 < -1.0)
            t0 = -1;
        end

        all0 = (lege.exev(t0,cfs));
        r0 = all0(1:dim);
        d0 = all0(dim+1:2*dim);
        d20 = all0(2*dim+1:end);

        rdiff0 = r0(:)-ref(:);
        dprime = rdiff0.'*d0(:);
        dprime2 = d0(:).'*d0(:) + d0(:).'*d20(:);

        if abs(dprime) < thresh
            newtsuccess = true;
            break;
        end

        if t0 == 1 && dprime < 0
            newtsuccess = true;
            break;
        end
        if t0 == -1 && dprime > 0
            newtsuccess = true;
            break;
        end

    end

    dist2newt = sum(abs(rdiff0.^2));

    if dist2newt > dist2all(i)
        newtsuccess = false;
    else
        ts(i) = t0;
        rs(:,i) = r0;
        ds(:,i) = d0;
        d2s(:,i) = d20;
        dist2s(i) = dist2newt;
    end

    if ~newtsuccess
        % try levenberg if Newton fails
        t0 = t(ipt(i));

        all0 = lege.exev(t0,cfs);
        r0 = all0(1:dim);
        d0 = all0(dim+1:2*dim);
        d20 = all0(2*dim+1:end);

        rdiff0 = r0(:) - ref;
        dprime = rdiff0.'*d0(:);
        dprime2 = d0(:).'*d0(:);
        lam = dprime2;

        dist0 = sum(abs(rdiff0).^2);

        for iii = 1:maxnewt

        % Levenberg step
        
            deltat = - dprime/(dprime2+lam);
    
            t1 = t0+deltat;
            
            if (t1 > 1.0)
                t1 = 1;
            end
            if (t1 < -1.0)
                t1 = -1;
            end
    
            all1 = (lege.exev(t1,cfs));
            r1 = all1(1:dim);
            d1 = all1(dim+1:2*dim);
            d21 = all1(2*dim+1:end);
    
            rdiff1 = r1(:)-ref(:);
            dist1 = sum(abs(rdiff1).^2);

            if dist1 > dist0
                lam = lam*2;
            else
                t0 = t1;
                r0 = r1;
                d0 = d1;
                d20 = d21;
                rdiff0 = rdiff1;
                dprime = rdiff0.'*d0(:);
                dprime2 = d0(:).'*d0(:) + d0(:).'*d20(:);
                dist0 = dist1;
                lam = lam/3;
                if abs(dprime) < thresh
                    break;
                end
        
                if t0 == 1 && dprime < 0
                    break;
                end
                if t0 == -1 && dprime > 0
                    break;
                end
            end
        end

        ts(i) = t0;
        rs(:,i) = r0;
        ds(:,i) = d0;
        d2s(:,i) = d20;
        dist2s(i) = dist0;
    end
        
end

