function [spmat] = buildmattd_fast(chnkr,kern,opdims,type,wts,ilist,logquad)
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
n = chnkr.n;
h = chnkr.h;
data = [];
if(chnkr.hasdata)
    data = chnkr.data;
end

xs1 = logquad.xs1;
wts1 = logquad.wts1;
xs0 = logquad.xs0;
wts0 = logquad.wts0;
ainterp1 = logquad.ainterp1;
ainterp1kron = logquad.ainterp1kron;
ainterps0 = logquad.ainterps0;
ainterps0kron = logquad.ainterps0kron;

mmat = k*nch*opdims(1); nmat = k*nch*opdims(2);
spmat = spalloc(mmat,nmat,k*nch*opdims(1)*k*3*opdims(2));

% overwrite nbor and self
for j = 1:nch

    jmat = 1 + (j-1)*k*opdims(2);
    jmatend = j*k*opdims(2);
    
    ibefore = adj(1,j);
    iafter = adj(2,j);

    % neighbors
    
    if ibefore > 0
      if ~isempty(ilist) && ismember(ibefore,ilist) && ismember(j,ilist) 
        % skip construction if both chunks are in the "bad" chunk list
      else  
        
        submat = chnk.quadggq.nearbuildmat(r,d,n,d2,h,data,ibefore,j, ...
           kern,opdims,xs1,wts1,ainterp1kron,ainterp1);    

        imat = 1 + (ibefore-1)*k*opdims(1);
        imatend = ibefore*k*opdims(1);

        spmat(imat:imatend,jmat:jmatend) = submat;
      end
    end
    
    if iafter > 0
      if ~isempty(ilist) && ismember(iafter,ilist) && ismember(j,ilist) 
        % skip construction if both chunks are in the "bad" chunk list
      else  
        
        submat = chnk.quadggq.nearbuildmat(r,d,n,d2,h,data,iafter,j, ...
           kern,opdims,xs1,wts1,ainterp1kron,ainterp1);
        

        imat = 1 + (iafter-1)*k*opdims(1);
        imatend = iafter*k*opdims(1);

        spmat(imat:imatend,jmat:jmatend) = submat;
      end
    end
    % self
    if ~isempty(ilist) && ismember(j,ilist) 
      % skip construction if the chunk is in the "bad" chunk list  
    else
      
      submat = chnk.quadggq.diagbuildmat(r,d,n,d2,h,data,j,kern,opdims,...
        xs0,wts0,ainterps0kron,ainterps0);
      

      imat = 1 + (j-1)*k*opdims(1);
      imatend = j*k*opdims(1);

      spmat(imat:imatend,jmat:jmatend) = submat;
    end
    
end
	 

end
