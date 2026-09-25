function [sysmat,auxquads] = buildmatparts(chnkr,kern,opdims,ilist,nonsmoothonly,corrections,parts,auxquads)
%CHNK.QUADGGQ.BUILDMATPARTS build matrix for given kernel and chnkr
% description of boundary, using special quadrature for self and neighbor
% panels. A different special quadrature is used for each element of parts.
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
%           kern is only used for the smooth (far) interactions.
%   opdims - (2) dimension of the kernel, for scalar kernels opdims(1:2) = 1;
%   ilist - cell array of integer arrays ([]), list of panel interactions that 
%          should be ignored when constructing matrix entries or quadrature
%          corrections. 
%   nonsmoothonly - boolean, if true only the self and neighbor panels are
%          computed and returned in a sparse matrix
%   corrections - boolean, if true (and nonsmoothonly) return the
%          corrections to the smooth quadrature rule instead
%   parts - struct with one kernel per singularity type, as in kernel.parts,
%          e.g. parts.log and parts.pv
%   auxquads - struct ([]), structure containing auxilliary quadrature
%          nodes, weights and related interpolation matrices for each
%          rule, e.g. auxquads.ggqlog (see chnk.quadggq.setup). Rules
%          that are missing or for a different order are recomputed.
%
% Ouput:
%   sysmat - the system matrix for discretizing integral operator whose kernel 
%            is defined by kern with a density on the domain defined by chnkr
%   auxquads - the input auxquads with any rules that were computed

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

% one kernel per quadrature rule, parts sharing a rule are added
temp = eye(opdims(2));
rules = struct();
for sname = fieldnames(parts)'
    type = sname{1};
    if strcmp(type, 'smooth')
        %TODO: make a reasonable method for smooth with removable
        type = 'log';
    end
    % recompute if missing or stored for a different order k
    if ~isfield(auxquads, ['ggq' type]) || size(auxquads.(['ggq' type]).ainterp1, 2) ~= k
        auxquads.(['ggq' type]) = chnk.quadggq.setup(k,type);
    end
    aq = auxquads.(['ggq' type]);

    parteval = parts.(sname{1}).eval;
    if isfield(rules, type)
        preveval = rules.(type).eval;
        parteval = @(s,t) preveval(s,t) + parteval(s,t);
    end
    rules.(type).eval = parteval;
    rules.(type).auxquads = aq;
    rules.(type).ainterp1kron = kron(aq.ainterp1,temp);
    rules.(type).ainterps0kron = cellfun(@(a) kron(a,temp), aq.ainterps0, ...
        'UniformOutput', false);
end
rulenames = fieldnames(rules)';

if corrections
    wtss = chnkr.wts;
    wtss = repmat(wtss(:).',opdims(2),1); wtss = reshape(wtss,opdims(2)*k,nch);
    indd = kron(eye(k),true(opdims(1),opdims(2)));
    indd = indd(:) > 0;
else
    wtss = [];
    indd = [];
end

mmat = k*nch*opdims(1); nmat = k*nch*opdims(2);

nnz = k*nch*opdims(1)*k*3*opdims(2);
nnz1 = k*opdims(1)*k*opdims(2);
v = zeros(nnz,1);
iind = zeros(nnz,1);
jind = zeros(nnz,1);
[jj1,ii1] = meshgrid(0:k*opdims(2)-1,0:k*opdims(1)-1);

% nbor and self, summing over the rules
ict = 0;
for j = 1:nch

    jmat = 1 + (j-1)*k*opdims(2);

    % neighbor before, neighbor after, self
    blocks = [adj(1,j), adj(2,j), j];
    for iblock = 1:3
        i = blocks(iblock);
        if i <= 0
            continue
        end
        % skip construction if both chunks are in the "bad" chunk list
        if ~isempty(ilist) && ismember(i,ilist) && ismember(j,ilist)
            continue
        end

        submat = 0;
        for rulename = rulenames
            p = rules.(rulename{1});
            if iblock == 3
                submat = submat + chnk.quadggq.diagbuildmat(r,d,n,d2,data,j,p.eval,opdims, ...
                    p.auxquads.xs0,p.auxquads.wts0,p.ainterps0kron,p.auxquads.ainterps0, ...
                    corrections,wtss,indd);
            else
                submat = submat + chnk.quadggq.nearbuildmat(r,d,n,d2,data,i,j,p.eval,opdims, ...
                    p.auxquads.xs1,p.auxquads.wts1,p.ainterp1kron,p.auxquads.ainterp1, ...
                    corrections,wtss);
            end
        end

        imat = 1 + (i-1)*k*opdims(1);
        iind(ict+(1:nnz1)) = ii1(:)+imat;
        jind(ict+(1:nnz1)) = jj1(:)+jmat;
        v(ict+(1:nnz1)) = submat(:);
        ict = ict + nnz1;
    end

end

iind = iind(1:ict);
jind = jind(1:ict);
v    = v(1:ict);

if nonsmoothonly
    sysmat = sparse(iind,jind,v,mmat,nmat);
else
    % do smooth weight for all, then overwrite nbor and self
    wts = chnkr.wstor;
    sysmat = chnk.quadnative.buildmat(chnkr,kern,opdims,1:nch,1:nch,wts);
    sysmat(sub2ind([mmat,nmat],iind,jind)) = v;
end

end
