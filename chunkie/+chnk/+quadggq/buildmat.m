function [sysmat] = buildmat(chnkr,kern,opdims,type)
%CHNK.QUADGGQ.BUILDMAT build matrix for given kernel and chnkr 
% description of boundary, using special quadrature for self
% and neighbor panels.
%
%  

k = chnkr.k;
nch = chnkr.nch;
r = chnkr.r;
adj = chnkr.adj;
d = chnkr.d;
d2 = chnkr.d2;
h = chnkr.h;

[~,wts,u] = lege.exps(k);

if strcmpi(type,'log')

    qavail = chnk.quadggq.logavail();
    [~,i] = min(abs(qavail-k));
    assert(qavail(i) == k,'order %d not found, consider using order %d chunks', ...
        k,qavail(i));
    [xs1,wts1,xs0,wts0] = chnk.quadggq.getlogquad(k);
else
    error('type not available')
end

ainterp1 = lege.matrin(k,xs1);
temp = eye(opdims(2));
ainterp1kron = kron(ainterp1,temp);

nquad0 = size(xs0,1);

ainterps0kron = zeros(opdims(2)*nquad0,opdims(2)*k,k);
ainterps0 = zeros(nquad0,k,k);

for j = 1:k
    xs0j = xs0(:,j);
    ainterp0_sm = lege.matrin(k,xs0j);
    ainterps0(:,:,j) = ainterp0_sm;
    ainterps0kron(:,:,j) = kron(ainterp0_sm,temp);
end

% do smooth weight for all
sysmat = chnk.quadnative.buildmat(chnkr,kern,opdims,1:nch,1:nch,wts);

% overwrite nbor and self
for j = 1:nch

    jmat = 1 + (j-1)*k*opdims(2);
    jmatend = j*k*opdims(2);
    
    ibefore = adj(1,j);
    iafter = adj(2,j);

    % neighbors
    
    if ibefore > 0
        submat = chnk.quadggq.nearbuildmat(r,d,d2,h,ibefore,j, ...
            kern,opdims,u,xs1,wts1,ainterp1kron,ainterp1);
    
        imat = 1 + (ibefore-1)*k*opdims(1);
        imatend = ibefore*k*opdims(1);

        sysmat(imat:imatend,jmat:jmatend) = submat;
    end
    
    if iafter > 0
        submat = chnk.quadggq.nearbuildmat(r,d,d2,h,iafter,j, ...
            kern,opdims,u,xs1,wts1,ainterp1kron,ainterp1);

        imat = 1 + (iafter-1)*k*opdims(1);
        imatend = iafter*k*opdims(1);
        
        sysmat(imat:imatend,jmat:jmatend) = submat;
    end
    
    % self
    
    submat = chnk.quadggq.diagbuildmat(r,d,d2,h,j,kern,opdims,...
        u,xs0,wts0,ainterps0kron,ainterps0);

    imat = 1 + (j-1)*k*opdims(1);
    imatend = j*k*opdims(1);

    sysmat(imat:imatend,jmat:jmatend) = submat;
    
end
	 

end
