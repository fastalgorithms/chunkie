function val = contract(comps,dim,vs)
%CHNK.RADPP.CONTRACT derivative of a radial function along given directions
%
% Syntax: val = chnk.radpp.contract(comps,dim,vs)
%
% comps is (nt,ns,ncomp), the derivative components with respect to the
% target from CHNK.RADPP.ROTATE, and vs is a cell array with one entry
% per derivative, each a cell array {vx,vy[,vz]} of direction components
% (scalars or arrays that expand to (nt,ns)). val is (nt,ns), the
% derivative tensor of order numel(vs) contracted with the directions,
%
%    val = sum_{i1...ik} D_{i1...ik} G v1_{i1} ... vk_{ik}.
%
% see also CHNK.RADPP.ROTATE, CHNK.RADPP.TUPLES

k = numel(vs);
if k == 0
    val = comps(:,:,1);
    return
end

% components of lower order come first
off = 0;
for j = 0:k-1
    off = off + nchoosek(dim+j-1,j);
end

tup = chnk.radpp.tuples(dim,k);
sz = [repmat(dim,1,k) 1];

val = 0;
for m = 1:dim^k
    sub = cell(1,k);
    [sub{:}] = ind2sub(sz,m);
    idx = [sub{:}];
    pos = find(all(tup == sort(idx),2));
    term = comps(:,:,off+pos);
    for j = 1:k
        term = term.*vs{j}{idx(j)};
    end
    val = val + term;
end

end
