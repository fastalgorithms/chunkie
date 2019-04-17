function [rc,dc,d2c] = chunkerexps(chnkr)
%CHUNKEREXPS compute Legendre coefficients for the position, derivative,
% and second derivative of the chunker object on each chunk.

k = chnkr.k; nch = chnkr.nch; dim = chnkr.dim;

[~,~,u] = legeexps(k);

rc = ipermute(reshape(u*reshape(permute(chnkr.r,[2 1 3]),k,nch*dim), ...
        k,dim,nch),[2 1 3]);

if nargout > 1
    dc = ipermute(reshape(u*reshape(permute(chnkr.d,[2 1 3]),k,nch*dim), ...
        k,dim,nch),[2 1 3]);
end

if nargout > 2 
    d2c = ipermute(reshape(u*reshape(permute(chnkr.d2,[2 1 3]),k,nch*dim), ...
        k,dim,nch),[2 1 3]);
end