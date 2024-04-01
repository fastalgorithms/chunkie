function [spmat] = buildmattd(chnkr,kern,opdims,type,auxquads,ilist)
%CHNK.QUADGGQ.BUILDMAT build matrix for given kernel and chnkr 
% description of boundary, using special quadrature for self
% and neighbor panels.
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
n = chnkr.n;

data = [];
if (chnkr.hasdata)
    data = chnkr.data;
end

if nargin < 4 || isempty(type)
    type = 'log';
end
if nargin < 5 || isempty(auxquads)
    auxquads = chnk.quadggq.setup(k,type);
end

temp = eye(opdims(2));

xs1 = auxquads.xs1;
wts1 = auxquads.wts1;
xs0 = auxquads.xs0;
wts0 = auxquads.wts0;

ainterp1 = auxquads.ainterp1;
ainterp1kron = kron(ainterp1,temp);

ainterps0 = auxquads.ainterps0;
ainterps0kron = cell(k,1);
for j = 1:k
    ainterps0kron{j} = kron(ainterps0{j},temp);
end

% do smooth weight for all
%sysmat = chnk.quadnative.buildmat(chnkr,kern,opdims,1:nch,1:nch,wts);
mmat = k*nch*opdims(1); nmat = k*nch*opdims(2);

nnz = k*nch*opdims(1)*k*3*opdims(2);
nz_found = 0;
nnz1 = k*opdims(1)*k*opdims(2);
v = zeros(nnz,1);
iind = zeros(nnz,1);
jind = zeros(nnz,1);

jmat = 0;
jmatend = k*opdims(2)-1;

imat = 0;
imatend = k*opdims(1)-1;
[jj1,ii1] = meshgrid(jmat:jmatend,imat:imatend);


% overwrite nbor and self
ict = 1;
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
            submat = chnk.quadggq.nearbuildmat(r,d,n,d2,data,ibefore,j, ...
                kern,opdims,xs1,wts1,ainterp1kron,ainterp1);
    
            imat = 1 + (ibefore-1)*k*opdims(1);
            induse = ict:ict+nnz1-1;
            iind(induse) = ii1(:)+imat;
            jind(induse) = jj1(:)+jmat;
            nz_found = nz_found + numel(induse);
            v(induse) = submat(:);
            ict = ict + nnz1;
        end
    end
    
    if iafter > 0
      if ~isempty(ilist) && ismember(iafter,ilist) && ismember(j,ilist) 
        % skip construction if both chunks are in the "bad" chunk list
      else      
        submat = chnk.quadggq.nearbuildmat(r,d,n,d2,data,iafter,j, ...
            kern,opdims,xs1,wts1,ainterp1kron,ainterp1);
        

        imat = 1 + (iafter-1)*k*opdims(1);
        induse = ict:ict+nnz1-1;
        iind(induse) = ii1(:)+imat;
        jind(induse) = jj1(:)+jmat;
        nz_found = nz_found + numel(induse);
        v(induse) = submat(:);
        ict = ict + nnz1;
      end
    end
    
    % self
    if ~isempty(ilist) && ismember(j,ilist) 
      % skip construction if the chunk is in the "bad" chunk list  
    else
      submat = chnk.quadggq.diagbuildmat(r,d,n,d2,data,j,kern,opdims,...
          xs0,wts0,ainterps0kron,ainterps0);

      imat = 1 + (j-1)*k*opdims(1);
      
      
      induse = ict:(ict+nnz1-1);
      iind(induse) = ii1(:)+imat;
      jind(induse) = jj1(:)+jmat;
      nz_found = nz_found + numel(induse);
      v(induse) = submat(:);
      ict = ict + nnz1;
    end
    
end

iind = iind(1:nz_found);
jind = jind(1:nz_found);
v    = v(1:nz_found);
spmat = sparse(iind,jind,v,mmat,nmat);
	 

end
