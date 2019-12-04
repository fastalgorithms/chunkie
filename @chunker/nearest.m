function [rn,dn,d2n,dist,tn,ichn] = nearest(chnkr,ref,ich,x,u)
%NEAREST Find nearest point on chunker for given reference point. 
%If ich provided, only checks chunks in ich. If not, checks all chunks.
%
%         k - number of points per chunk
%         xnodes - the k Legendre points
%         u - the matrix created by legewhts mapping points to coefs
%         chunk - points on the chunk at legendre nodes
%         der - derivatices w.r.t. "some" parameterization
%         der2 - 2nd derivatives w.r.t. the same parameterization
%         h - normalizing factor for the parameterization
%         targ - the distance from chunk ich to targ is desired

%       output:
%         dist - the distance from the chunk to targ
%         xy - the point lying on the chunk that is the closest, this
%             point is determined by newton
%         xy_norm - the outward normal, i.e. into the unbounded domain,
%             assuming that the chunk is parameterized counter-clockwise


%       find the closest legendre node



%        ifwhts = 0
%        call legewhts(k, xnodes, whts, ifwhts)

if nargin < 4
    [x,w,u] = lege.exps(chnkr.k);
end

dist = Inf;

if nargin < 3
    ich = 1:nch;
end

for i = 1:length(ich)
    % grab chunk data
    ii = ich(i);
    r = chnkr.r(:,:,ii);
    d = chnkr.d(:,:,ii);
    d2 = chnkr.d2(:,:,ii);
    h = chnkr.h(ii);
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
    maxnewt = 200;
    thresh = 1.0d-10;
    thresh = 1.0d-9;
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

