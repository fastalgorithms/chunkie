function [sysmat] = buildmat(chnkr,kern,opdims,type,auxquads,ilist)
%CHNK.QUADGGQ.BUILDMAT build matrix for given kernel and chnkr 
% description of boundary, using special quadrature for self
% and neighbor panels.
%
% Input:
%   chnkr - chunker object describing boundary
%   kern  - kernel function. By default, this should be a function handle
%           accepting input of the form kern(srcinfo,targinfo), where srcinfo
%           and targinfo are in the ptinfo struct format, i.e.
%                ptinfo.r - positions (2,:) array
%                ptinfo.d - first derivative in underlying
%                     parameterization (2,:)
%                ptinfo.n - unit normals (2,:)
%                ptinfo.d2 - second derivative in underlying
%                     parameterization (2,:)
%   opdims - (2) dimension of the kernel, for scalar kernels opdims(1:2) = 1;
%
% Optional input: quantities in brackets indicate default settings
%  type - string ('log'), type of singularity of kernel. Type
%          can take on the following arguments:
%             log => logarithmically singular kernels
%             pv => principal value singular kernels
%             hs => hypersingular kernels
%  auxquads - struct (chnk.quadggq.setuplogquads), structure containing
%             auxilliary quadrature nodes, weights and related
%             interpolation matrices.
%  ilist - cell array of integer arrays ([]), list of panel interactions that 
%          should be ignored when constructing matrix entries or quadrature
%          corrections. 
%
% Ouput:
%   sysmat - the system matrix for discretizing integral operator whose kernel 
%            is defined by kern with a density on the domain defined by chnkr
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
wts = chnkr.wstor;
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
            submat = chnk.quadggq.nearbuildmat(r,d,n,d2,data,ibefore,j, ...
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
        submat = chnk.quadggq.nearbuildmat(r,d,n,d2,data,iafter,j, ...
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
      submat = chnk.quadggq.diagbuildmat(r,d,n,d2,data,j,kern,opdims,...
          xs0,wts0,ainterps0kron,ainterps0);

      imat = 1 + (j-1)*k*opdims(1);
      imatend = j*k*opdims(1);

      sysmat(imat:imatend,jmat:jmatend) = submat;
    end
    
end
	 

end
