function rctr = ctr(obj)
%CTR center of chunkgraph object
%
% This is the average of cgrph.max and cgrph.min
%
% Syntax: rctr = ctr(cgrph)
%
% Input:
%   chnkr - chunkgraph object
%
% Output:
%   rctr - chnkr.dim length vector, average of max and min of chnkr
%
% Examples:
%   rctr = ctr(cgrph);
%

% author: Tristan Goodwill

rctr = (max(obj) + min(obj))/2;
