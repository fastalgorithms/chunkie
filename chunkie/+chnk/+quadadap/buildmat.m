function [sysmat] = buildmat(chnkr,kern,opdims,type,opts)
%CHNK.QUADADAP.BUILDMAT build matrix for given kernel and chnkr 
% description of boundary, using special quadrature for self panel
% and adaptive quadrature for neighbor panels. Optionally, adaptive 
% quadrature can be applied to nearly singular interactions, i.e. 
% targets within a chunk length of each chunk
%
%  

if nargin < 5
    opts = [];
end

robust = false;
if isfield(opts,'robust')
    robust = opts.robust;
end

k = chnkr.k;
nch = chnkr.nch;
r = chnkr.r;
adj = chnkr.adj;
d = chnkr.d;
d2 = chnkr.d2;
h = chnkr.h;
n = chnkr.n;

[t,wts,u] = lege.exps(k);
bw = lege.barywts(k);

k2 = max(27,k+1);
[t2,w2] = lege.exps(k2);

if strcmpi(type,'log')

    qavail = chnk.quadggq.logavail();
    [~,i] = min(abs(qavail-k));
    assert(qavail(i) == k,'order %d not found, consider using order %d chunks', ...
        k,qavail(i));
    [~,~,xs0,wts0] = chnk.quadggq.getlogquad(k);
else
    error('type not available')
end

nquad0 = size(xs0,1);

ainterps0kron = zeros(opdims(2)*nquad0,opdims(2)*k,k);
ainterps0 = zeros(nquad0,k,k);

temp = eye(opdims(2));

for i = 1:k
    xs0j = xs0(:,i);
    ainterp0_sm = lege.matrin(k,xs0j);
    ainterps0(:,:,i) = ainterp0_sm;
    ainterps0kron(:,:,i) = kron(ainterp0_sm,temp);
end

% do smooth weight for all
sysmat = chnk.quadnative.buildmat(chnkr,kern,opdims,1:nch,1:nch,wts);

% overwrite nbor and self
for i = 1:nch

    jmat = 1 + (i-1)*k*opdims(2);
    jmatend = i*k*opdims(2);
    
    ibefore = adj(1,i);
    iafter = adj(2,i);

    % neighbors
    
    if ibefore > 0
        rt = r(:,:,ibefore);
        dt = d(:,:,ibefore);
        d2t = d2(:,:,ibefore);
        nt = n(:,:,ibefore);
        submat = chnk.adapgausswts(r,d,n,d2,h,t,bw,i,rt,dt,nt,d2t, ...
            kern,opdims,t2,w2);

        imat = 1 + (ibefore-1)*k*opdims(1);
        imatend = ibefore*k*opdims(1);

        sysmat(imat:imatend,jmat:jmatend) = submat;
    end
    if iafter > 0
        rt = r(:,:,iafter);
        dt = d(:,:,iafter);
        nt = n(:,:,iafter);
        d2t = d2(:,:,iafter);
        submat = chnk.adapgausswts(r,d,n,d2,h,t,bw,i,rt,dt,nt,d2t, ...
            kern,opdims,t2,w2);

        imat = 1 + (iafter-1)*k*opdims(1);
        imatend = iafter*k*opdims(1);

        sysmat(imat:imatend,jmat:jmatend) = submat;
    end
    
    % self
    
    submat = chnk.quadggq.diagbuildmat(r,d,n,d2,h,[],i,kern,opdims,...
        xs0,wts0,ainterps0kron,ainterps0);

    imat = 1 + (i-1)*k*opdims(1);
    imatend = i*k*opdims(1);

    sysmat(imat:imatend,jmat:jmatend) = submat;
    
end
	
if robust
   
    % pair-wise distances of chunk pts
    
    dim = chnkr.dim;
    k = chnkr.k;
    nch = chnkr.nch;
    
    dists = zeros(k*nch,k*nch);
    for i = 1:dim
        ri = r(i,:);
        dists = dists + (ri - ri.').^2;
    end
    dists = sqrt(dists);
    
    wtschnk = weights(chnkr);
    chnklen = sum(wtschnk,1);

    % do adaptive quadrature for points that aren't
    % neighbors/self points but are too close
    
    for i = 1:nch
        istart = (i-1)*k+1;
        iend = i*k;
        jmat = 1 + (i-1)*k*opdims(2);
        jmatend = i*k*opdims(2);

        chnkdisti = dists(:,istart:iend);
        chnkdistmini = min(chnkdisti,[],2);
        targfix = find(chnkdistmini < chnklen(i));
        
        iright = adj(2,i);
        ileft = adj(1,i);
        targignore = istart:iend;
        if iright > 0
            ir1 = (iright-1)*k+1;
            ir2 = iright*k;
            targignore = [targignore, ir1:ir2];
        end
        if ileft > 0
            il1 = (ileft-1)*k+1;
            il2 = ileft*k;
            targignore = [targignore, il1:il2];
        end
        
        targfix = setdiff(targfix,targignore(:));
        
        rt = r(:,targfix);
        dt = d(:,targfix);
        d2t = d2(:,targfix);
        nt = n(:,targfix);

        submat = chnk.adapgausswts(r,d,n,d2,h,t,bw,i,rt,dt,nt,d2t, ...
                kern,opdims,t2,w2);
            
        imats = bsxfun(@plus,(1:opdims(1)).',opdims(1)*(targfix(:)-1).');
        imats = imats(:);
        sysmat(imats,jmat:jmatend) = submat;
            
    end
    
    
end

end
