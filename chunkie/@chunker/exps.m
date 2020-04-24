function [rc,dc,d2c] = exps(chnkr)
%EXPS compute Legendre coefficients for the position, derivative,
% and second derivative of the chunker object on each chunk.
%
% Syntax: [rc,dc,d2c] = exps(chnkr)
%
% Input:
%   chnkr - chunker object
%
% Output:
%   rc - coefficients for position on each chunk
%   dc - coefficients for derivative of position on each chunk
%   d2c - coefficients for 2nd derivative of position on each chunk
%
% Examples:
%   [rc,dc,d2c] = exps(chnkr)
%   [~,dc] = exps(chnkr)
%

% author: Travis Askham (askhamwhat@gmail.com)

k = chnkr.k; nch = chnkr.nch; dim = chnkr.dim;

[~,~,u] = lege.exps(k);

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
