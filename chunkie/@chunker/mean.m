function rmean = mean(obj)
%MEAN weighted average position over chunker nodes
%
% This is the average of chnkr.r weighted by the integration weights,
% i.e. an approximation to the average position of the curve
%
% Syntax: rmean = mean(chnkr)
%
% Input:
%   chnkr - chunker object
%
% Output:
%   rmean - chnkr.dim length vector, weighted average of node positions
%
% Examples:
%   rmean = mean(chnkr);
%

% author: Tristan Goodwill

rs = reshape(real(obj.r),obj.dim,obj.k*obj.nch);
wts = obj.wts(:).';

rmean = sum(rs.*wts,2)/sum(wts);
