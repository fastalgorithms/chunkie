function errs = chunk_fun_error(cg, fval)
%CHUNK_FUN_ERROR estimate resolution of a function on each chunk of a chunkgraph
%
% Applies chunker/chunk_fun_error to each edge and concatenates results.
% The function values fval should be ordered consistently with cg.r, i.e.
% edge by edge in the order of cg.echnks.
%
% Syntax: errs = chunk_fun_error(cg, fval)
%
% Input:
%   cg   - chunkgraph object
%   fval - nfuns x npt, nfuns x k x nch, or vector of length nfuns*npt
%          array of function values at the chunkgraph nodes (cg.r ordering)
%
% Output:
%   errs - nfuns x nch array; errs(j,i) is the estimated error of
%          function j on chunk i (max abs value of last two Legendre
%          series coefficients), with chunks ordered edge by edge
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

errs  = zeros(nfuns, nch);
ipt   = 0;
ich   = 0;
for j = 1:length(cg.echnks)
    nchj  = cg.echnks(j).nch;
    nptj  = k * nchj;
    fvalj = fval(:, ipt+1:ipt+nptj);
    errs(:, ich+1:ich+nchj) = chunk_fun_error(cg.echnks(j), fvalj);
    ipt = ipt + nptj;
    ich = ich + nchj;
end

end
