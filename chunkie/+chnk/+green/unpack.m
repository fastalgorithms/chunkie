function c = unpack(v,dim,k)
%CHNK.GREEN.UNPACK derivative components of order k from a Green's
% function output
%
% Syntax: c = chnk.green.unpack(v,dim,k)
%
% v is the k-th derivative output of a Green's function, either packed
% (nt,ns,ncomp), with the components in lexicographic order, or
% full (nt,ns,dim,...,dim) with k trailing dimensions. c is (nt*ns,ncomp),
% the packed form.
%
% see also KERNEL.GREEN, CHNK.GREEN.KERN

tup = nchoosek(1:dim+k-1,k) - (0:k-1);
ncomp = size(tup,1);

if ndims(v) <= 3
    if size(v,3) ~= ncomp
        error('CHNK.GREEN.UNPACK: derivative %d should have %d components.', ...
            k,ncomp);
    end
    c = reshape(v,[],ncomp);
    return
end

sz = size(v);
if ndims(v) ~= 2+k || any(sz(3:2+k) ~= dim)
    error('CHNK.GREEN.UNPACK: derivative %d has the wrong size.',k);
end

v = reshape(v,size(v,1)*size(v,2),dim^k);
c = zeros(size(v,1),ncomp);
for m = 1:ncomp
    sub = num2cell(tup(m,:));
    c(:,m) = v(:,sub2ind([repmat(dim,1,k) 1],sub{:}));
end

end
