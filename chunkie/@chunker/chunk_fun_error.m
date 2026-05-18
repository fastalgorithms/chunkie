function errs = chunk_fun_error(chnkr, fval)
%CHUNK_FUN_ERROR estimate resolution of a function on each chunk
%
% For each chunk, returns the maximum absolute value of the last two
% Legendre coefficients as a proxy for the truncation error.
%
% Syntax: errs = chunk_fun_error(chnkr, fval)
%
% Input:
%   chnkr  - chunker object
%   fval   - array of function values at the chunker nodes with nfuns *
%            npts elements
%
% Output:
%   errs   - nfuns x nch array; errs(j,i) is the estimated error of
%            function j on chunk i
%
% see also CHUNKER.REFINE

% author: Tristan Goodwill

k    = chnkr.k;
nch  = chnkr.nch;

% accept:
%   - vector of length nfuns*k*nch (any orientation)
%   - nfuns x k*nch matrix
%   - nfuns x k x nch array
npts = k*nch;
if isvector(fval)
    nfuns = numel(fval) / npts;
    assert(nfuns == round(nfuns), ...
        'CHUNKER:chunk_fun_error length of fval must be divisible by k*nch');
    fval = reshape(fval, nfuns, npts);
else
    nfuns = size(fval, 1);
end

% u maps k function values at Legendre nodes to k Legendre coefficients
[~,~,u] = lege.exps(k);

% accept nfuns x k*nch or nfuns x k x nch
fval = reshape(fval, nfuns, k, nch);
errs = zeros(nfuns, nch);
for i = 1:nch
    % coefs is k x nfuns
    coefs = u * fval(:,:,i).';
    % last two rows give the two highest-degree coefficients
    errs(:,i) = max(abs(coefs(end-1:end, :)), [], 1);
end

end
