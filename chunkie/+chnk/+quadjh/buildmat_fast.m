function [sysmat] = buildmat_fast(chnkr,kern,opdims,type,wts,ilist,logquad)
%CHNK.QUADJH.BUILDMAT build matrix for given kernel and chnkr 
% description of boundary, using special quadrature for self
% and neighbor panels.
%
%  difference between this function and the function buildmat:
% all quadrature nodes, weights, and interpolation matrices for
% logarithmically and nearly logarithmically singular integrals 
% are precomputed and supplied as inputs.

k = chnkr.k;
nch = chnkr.nch;
r = chnkr.r;
adj = chnkr.adj;
d = chnkr.d;
n = chnkr.n;
d2 = chnkr.d2;
h = chnkr.h;

xs1 = logquad.xs1;
wts1 = logquad.wts1;
ainterp1 = logquad.ainterp1;
ainterp1kron = logquad.ainterp1kron;

%[~,wts] = lege.exps(k);
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
      if ~isempty(ilist) && ismember(ibefore,ilist) && ismember(j,ilist) 
        % skip correction if both chunks are in the "bad" chunk list        
      else
        submat = chnk.quadggq.nearbuildmat(r,d,n,d2,h,ibefore,j, ...
            kern,opdims,xs1,wts1,ainterp1kron,ainterp1);
    
        imat = 1 + (ibefore-1)*k*opdims(1);
        imatend = ibefore*k*opdims(1);

        sysmat(imat:imatend,jmat:jmatend) = submat;
      end
    end
    
    if iafter > 0
      if ~isempty(ilist) && ismember(iafter,ilist) && ismember(j,ilist) 
        % skip correction if both chunks are in the "bad" chunk list        
      else      
        submat = chnk.quadggq.nearbuildmat(r,d,n,d2,h,iafter,j, ...
            kern,opdims,xs1,wts1,ainterp1kron,ainterp1);

        imat = 1 + (iafter-1)*k*opdims(1);
        imatend = iafter*k*opdims(1);
        
        sysmat(imat:imatend,jmat:jmatend) = submat;
      end
    end
    
    % self
    if ~isempty(ilist) && ismember(j,ilist) 
      % skip correction if the chunk is in the "bad" chunk list        
    else
      submat = chnk.quadjh.diagbuildmat(r,d,n,d2,h,j,kern,opdims,wts,logquad);

      imat = 1 + (j-1)*k*opdims(1);
      imatend = j*k*opdims(1);

      sysmat(imat:imatend,jmat:jmatend) = submat;
    end
end
	 

end