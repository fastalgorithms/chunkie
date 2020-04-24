function lens = chunklen(chnkr,ich,w)
%CHUNKLEN lengths of chunks in chunker object
%
% Syntax: lens = chunklen(chnkr,ich,w)
%
% Input:
%   chnkr - chunker object
%
% Optional input:
%   ich - subset of chunks to get lengths of 
%   w - precomputed Legendre weights of order chnkr.k
%
% Output:
%   lens - length of chunks
%
% Examples:
%   lens = chunklen(chnkr);
%   ich = [1,3,9]; [~,w] = lege.exps(chnkr.k);
%   lens = chunklen(chnkr,ich,w);
%

% author: Travis Askham (askhamwhat@gmail.com)

if nargin < 2
    ich = 1:chnkr.nch;
end

k = chnkr.k;

if nargin < 3
  [~,w] = lege.exps(k);
end

nch1 = numel(ich);
wts = reshape(sqrt(sum((chnkr.d(:,:,ich(:))).^2,1)),k,nch1);
wts = wts.*bsxfun(@times,w(:),(chnkr.h(ich(:))).');
lens = sum(wts,1); lens = lens(:);

end
