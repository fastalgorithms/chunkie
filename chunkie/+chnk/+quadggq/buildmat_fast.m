function [sysmat] = buildmat_fast(chnkr,kern,opdims,type,wts,ilist,...
  xs1,wts1,xs0,wts0,ainterp1,ainterp1kron,ainterps0,ainterps0kron)
%CHNK.QUADGGQ.BUILDMAT build matrix for given kernel and chnkr 
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
        
    else
      submat = chnk.quadggq.diagbuildmat(r,d,n,d2,h,j,kern,opdims,...
          xs0,wts0,ainterps0kron,ainterps0);

      imat = 1 + (j-1)*k*opdims(1);
      imatend = j*k*opdims(1);

      sysmat(imat:imatend,jmat:jmatend) = submat;
    end
end
	 

end
