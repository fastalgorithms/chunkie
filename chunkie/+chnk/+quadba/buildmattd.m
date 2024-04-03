function [spmat] = buildmattd(chnkr,kern,quadorder,opdims,type)
%BUILDMATTD build sparse matrix corresponding to self and neighbor
% interactions for given kernel and chnkr description of boundary
%
%  


k = chnkr.k;
nch = chnkr.nch;
r = chnkr.r;
adj = chnkr.adj;
d = chnkr.d;
d2 = chnkr.d2;

[~,~,u] = lege.exps(k);

if strcmpi(type,'log')

    qavail = chnk.quadggq.logquadavail();
    [~,i] = min(abs(qavail-quadorder));
    if (qavail(i) ~= quadorder)
        warning('order %d not found, using order %d', ...
            quadorder,qavail(i));
        quadorder = qavail(i);
    end    
    [xs1,whts1,xs0,whts0] = chnk.quadggq.getlogquad(quadorder);
else
    error('type not available')
end

ainterp1_sm = lege.matrin(k,xs1);
temp = eye(opdims(2));
ainterp1 = kron(ainterp1_sm,temp);

nquad0 = size(xs0,1);

ainterps0 = zeros(opdims(2)*nquad0,opdims(2)*k,k);

for j = 1:k
    xs0j = xs0(:,j);
    ainterp0_sm = lege.matrin(k,xs0j);
    ainterps0(:,:,j) = kron(ainterp0_sm,temp);
end

mmat = k*nch*opdims(1); nmat = k*nch*opdims(2);
spmat = spalloc(mmat,nmat,k*nch*opdims(1)*k*3*opdims(2));

% overwrite nbor and self
for j = 1:nch

    jmat = 1 + (j-1)*k*opdims(2);
    jmatend = j*k*opdims(2);
    
    ibefore = adj(1,j);
    iafter = adj(2,j);

    % neighbors
    
    submat = chnk.quadggq.nearbuildmat(r,d,h,ibefore,j, ...
        kern,opdims,u,xs1,whts1,ainterp1);
    
    submat(submat == 0.0) = 1e-300; %hack
    
    imat = 1 + (ibefore-1)*k*opdims(1);
    imatend = ibefore*k*opdims(1);

    spmat(imat:imatend,jmat:jmatend) = submat;
    
    submat = chnk.quadggq.nearbuildmat(r,d,h,iafter,j, ...
        kern,opdims,u,xs1,whts1,ainterp1);
    
    submat(submat == 0.0) = 1e-300;
    
    imat = 1 + (iafter-1)*k*opdims(1);
    imatend = iafter*k*opdims(1);

    spmat(imat:imatend,jmat:jmatend) = submat;
    
    % self
    
    submat = chnk.quadggq.diagbuildmat(r,d,h,j,kern,opdims,...
        u,xs0,whts0,ainterps0);
    
    submat(submat == 0.0) = 1e-300;

    imat = 1 + (j-1)*k*opdims(1);
    imatend = j*k*opdims(1);

    spmat(imat:imatend,jmat:jmatend) = submat;
    
end
	 

end
