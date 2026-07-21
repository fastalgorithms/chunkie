function rctr = ctr(obj)
%CTR center of chunker object
%
% This is the average of chnkr.max and chnkr.min
%
% Syntax: rctr = ctr(chnkr)
%
% Input:
%   chnkr - chunker object
%
% Output:
%   rctr - chnkr.dim length vector, average of max and min of chnkr
%
% Examples:
%   rctr = ctr(chnkr);
%

% author: Tristan Goodwill

rctr = (max(obj) + min(obj))/2;
