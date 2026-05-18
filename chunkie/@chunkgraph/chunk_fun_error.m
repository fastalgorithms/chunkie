function errs = chunk_fun_error(cg, fval)
%CHUNK_FUN_ERROR estimate resolution of a function on each chunk of a chunkgraph
%
% For each chunk, returns the maximum absolute value of the last two
% Legendre coefficients as a proxy for the truncation error.
%
% Syntax: errs = chunk_fun_error(cg, fval)
%
% Input:
%   cg   - chunkgraph object
%   fval   - array of function values at the chunkgraph nodes with nfuns *
%            npts elements
%
% Output:
%   errs - nfuns x nch array; errs(j,i) is the estimated error of
%          function j on chunk i
%
% see also CHUNKER/CHUNK_FUN_ERROR, CHUNKGRAPH/REFINE

% author: Tristan Goodwill

k   = cg.k;
nch = sum([cg.echnks.nch]);

% handle vector input: infer nfuns from total length
if isvector(fval)
    nfuns = numel(fval) / (k * nch);
    assert(nfuns == round(nfuns), ...
        'CHUNKGRAPH:chunk_fun_error length of fval must be divisible by k*nch');
    fval = reshape(fval, nfuns, k * nch);
else
    nfuns = size(fval, 1);
    fval  = reshape(fval, nfuns, k * nch);
end

% u maps k function values at Legendre nodes to k Legendre coefficients
[~,~,u] = lege.exps(k);

fval = reshape(fval, nfuns, k, nch);
errs = zeros(nfuns, nch);
for i = 1:nch
    coefs = u * fval(:,:,i).';
    errs(:,i) = max(abs(coefs(end-1:end, :)), [], 1);
end

end
