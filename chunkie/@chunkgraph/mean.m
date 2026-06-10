function rmean = mean(obj)
%MEAN weighted average position over chunkgraph nodes
%
% This is the average of cgrph.r weighted by the integration weights,
% i.e. an approximation to the average position of the curve
%
% Syntax: rmean = mean(cgrph)
%
% Input:
%   chnkr - chunkgraph object
%
% Output:
%   rmean - cgrph.dim length vector, weighted average of node positions
%
% Examples:
%   rmean = mean(cgrph);
%

% author: Tristan Goodwill

rs = reshape(real(obj.r),obj.dim,obj.npt);
wts = obj.wts(:).';

rmean = sum(rs.*wts,2)/sum(wts);
