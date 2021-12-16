function [sysmat] = buildmat(chnkr,kern,opdims,type,auxquads,ilist)
%CHNK.QUADGGQ.BUILDMAT build matrix for given kernel and chnkr 
% description of boundary, using special quadrature for self
% and neighbor panels.
%
%  

% type: the type of singularity present in the kernel
% 
% NB: if type and auxquads are both present then auxquads will be used 
%       independent of type.

if (nargin < 3)
    error('not enough arguments in chnk.quadggq.buildmat');
end

if (nargin <6)
    ilist = [];
end

k = chnkr.k;
nch = chnkr.nch;
r = chnkr.r;
adj = chnkr.adj;
d = chnkr.d;
d2 = chnkr.d2;
h = chnkr.h;
n = chnkr.n;

data = [];
if (chnkr.hasdata)
    data = chnkr.data;
end

[~,wts] = lege.exps(k);

if (nargin == 4)
    if strcmpi(type,'log')
        auxquads = chnk.quadggq.setuplogquad(k,opdims);
    end
end 

if (nargin<4)
     auxquads = chnk.quadggq.setuplogquad(k,opdims);
end

%if strcmpi(type,'log')

    %qavail = chnk.quadggq.logavail();
    %[~,i] = min(abs(qavail-k));
    %assert(qavail(i) == k,'order %d not found, consider using order %d chunks', ...
    %    k,qavail(i));
    %[xs1,wts1,xs0,wts0] = chnk.quadggq.getlogquad(k);
    
    xs1 = auxquads.xs1;
    wts1 = auxquads.wts1;
    xs0 = auxquads.xs0;
    wts0 = auxquads.wts0;

    
    
%%%else
%%%%    error('type not available')
%%%end

%ainterp1 = lege.matrin(k,xs1);
ainterp1 = auxquads.ainterp1;

%temp = eye(opdims(2));
%ainterp1kron = kron(ainterp1,temp);
ainterp1kron = auxquads.ainterp1kron;

%nquad0 = size(xs0,1);

%%%%%ainterps0kron = zeros(opdims(2)*nquad0,opdims(2)*k,k);
%%%%%ainterps0 = zeros(nquad0,k,k);
%%%%%ainterps0 = logquad.ainterps0;
%%%%ainterps0kron = logquad.ainterps0kron;

%%%%%for j = 1:k
%%%%%    xs0j = xs0(:,j);
%%%%%    ainterp0_sm = lege.matrin(k,xs0j);
%%%%%    ainterps0(:,:,j) = ainterp0_sm;
%%%%%    ainterps0kron(:,:,j) = kron(ainterp0_sm,temp);
%%%%%end

ainterps0 = auxquads.ainterps0;
ainterps0kron = auxquads.ainterps0kron;

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
        % skip construction if both chunks are in the "bad" chunk list
        else
            submat = chnk.quadggq.nearbuildmat(r,d,n,d2,h,data,ibefore,j, ...
                kern,opdims,xs1,wts1,ainterp1kron,ainterp1);
    
            imat = 1 + (ibefore-1)*k*opdims(1);
            imatend = ibefore*k*opdims(1);

            sysmat(imat:imatend,jmat:jmatend) = submat;
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
        
        sysmat(imat:imatend,jmat:jmatend) = submat;
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

      sysmat(imat:imatend,jmat:jmatend) = submat;
    end
    
end
	 

end
