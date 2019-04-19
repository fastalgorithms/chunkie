function [sysmat] = chunkskernmat(chnkr,fkern,opdims,intparams)
%CHUNKSKERNMAT build matrix for given kernel and chnkr description of 
% boundary
%
%  

k = chnkr.k;
nch = chnkr.nch;
r = chnkr.r;
adj = chnkr.adj;
d = chnkr.d;
d2 = chnkr.d2;
h = chnkr.h;

dim = chnkr.dim;

[~,whts,u] = lege.exps(k);
intorder = intparams.intorder;
[xs1,whts1,xs0,whts0] = quad.brem.getquad(intorder);

%sysmat = zeros(opdims(1)*k*nch,opdims(2)*k*nch);

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

% do smooth weight for all
sysmat = chunksfarbuildmat(r,d,h,1:nch,1:nch,fkern,opdims,whts);

% overwrite nbor and self
for j = 1:nch

    jmat = 1 + (j-1)*k*opdims(2);
    jmatend = j*k*opdims(2);
    
    ibefore = adj(1,j);
    iafter = adj(2,j);

    % neighbors
    
    submat = chunksnearbuildmat(r,d,h,ibefore,j, ...
        fkern,opdims,u,xs1,whts1,ainterp1);
    
    imat = 1 + (ibefore-1)*k*opdims(1);
    imatend = ibefore*k*opdims(1);

    sysmat(imat:imatend,jmat:jmatend) = submat;
    
    submat = chunksnearbuildmat(r,d,h,iafter,j, ...
        fkern,opdims,u,xs1,whts1,ainterp1);
    
    imat = 1 + (iafter-1)*k*opdims(1);
    imatend = iafter*k*opdims(1);

    sysmat(imat:imatend,jmat:jmatend) = submat;
    
    % self
    
    submat = chunksdiagbuildmat(r,d,h,j,fkern,opdims,...
        u,xs0,whts0,ainterps0);

    imat = 1 + (j-1)*k*opdims(1);
    imatend = j*k*opdims(1);

    sysmat(imat:imatend,jmat:jmatend) = submat;
    
%     
%     jmat = 1 + (j-1)*k*opdims(2);
%     jmatend = j*k*opdims(2);
%     for i = 1:nch
%         imat = 1 + (i-1)*k*opdims(1);
%         imatend = i*k*opdims(1);
%         ibefore = adj(1,i);
%         iafter = adj(2,i);
% 
%         if (j == ibefore || j == iafter)
%             submat = chunksnearbuildmat(r,d,h,i,j, ...
%                 fkern,opdims,u,xs1,whts1,ainterp1);
%         elseif (j == i)
%             submat = chunksdiagbuildmat(r,d,h,j,fkern,opdims,...
%                 u,xs0,whts0,ainterps0);
%         else
%             submat = chunksfarbuildmat(r,d,h,i,j,fkern,opdims,whts);
%         end
%         
%         sysmat(imat:imatend,jmat:jmatend) = submat;
%     end
end
	 

end
