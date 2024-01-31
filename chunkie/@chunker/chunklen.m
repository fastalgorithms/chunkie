function lens = chunklen(chnkr,ich)
%CHUNKLEN lengths of chunks in chunker object
%
% Syntax: lens = chunklen(chnkr,ich,w)
%
% Input:
%   chnkr - chunker object
%
% Optional input:
%   ich - subset of chunks to get lengths of 
%
% Output:
%   lens - length of chunks
%
% Examples:
%   lens = chunklen(chnkr);
%   ich = [1,3,9]; 
%   lens = chunklen(chnkr,ich);
%

% author: Travis Askham (askhamwhat@gmail.com)

if nargin < 2
    ich = 1:chnkr.nch;
end

wts = chnkr.wts(:,ich);
lens = sum(wts,1); lens = lens(:);

end
