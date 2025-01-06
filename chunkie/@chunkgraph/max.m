function rmax = max(obj)
%MAX maximum value of each coordinate over chunkgraph nodes
% 
% This is *not* the maximum extent of the curve
%
% Syntax: rmax = max(chnkr)
%
% Input:
%   chnkr - chunker object
%
% Output:
%   rmax - chnkr.dim length vector of maximum node value in each direction
%
% Examples:
%   rmax = max(chnkr);
%

% author: Travis Askham (askhamwhat@gmail.com)

rmax = max(real(reshape(obj.r,obj.dim,obj.npt)),[],2);