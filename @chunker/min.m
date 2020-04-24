function rmin = min(obj)
%MIN minimum value of each coordinate over chunker nodes
% 
% This is *not* the minimum extent of the curve
%
% Syntax: rmin = min(chnkr)
%
% Input:
%   chnkr - chunker object
%
% Output:
%   rmin - chnkr.dim length vector of minimum node value in each direction
%
% Examples:
%   rmin = min(chnkr);
%

% author: Travis Askham (askhamwhat@gmail.com)

rmin = min(reshape(obj.r,obj.dim,obj.k*obj.nch),[],2);