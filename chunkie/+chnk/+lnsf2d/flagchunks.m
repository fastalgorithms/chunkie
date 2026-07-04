function flags = flagchunks(chnkr, prm, tau)
%CHNK.LNSF2D.FLAGCHUNKS flag chunks that are long compared to the
% boundary-layer decay length of the thermoviscous kernels, for which the
% standard GGQ self/neighbor rules lose accuracy (the smooth factors of
% the kernel vary on the scale 1/|k_h|, 1/|k_v| and are no longer
% resolved by degree k-1 polynomials on the chunk).
%
% Syntax: flags = chnk.lnsf2d.flagchunks(chnkr, prm, tau)
%
% Input:
%   chnkr - chunker object
%   prm - parameter struct from chnk.lnsf2d.params
%   tau - threshold (default 3.0): chunk j is flagged if
%         min(Im kh, Im kv) * chunklen(j) > tau
%
% Output:
%   flags - row vector of flagged chunk indices
%
% see also CHNK.LNSF2D.NEARCORRECT

if nargin < 3, tau = 3.0; end
kim = min(imag(prm.kh), imag(prm.kv));
lens = chunklen(chnkr);
flags = find(lens(:).'*kim > tau);

end
