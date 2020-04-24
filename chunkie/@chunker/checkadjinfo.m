function ier = checkadjinfo(chnkr)
%CHECKADJINFO checks adjacency vector for errors. based on attempting
% a sort.
%
% This is merely a wrapper for sortinfo.
%
% Syntax: ier = checkadjinfo(chnkr)
%
% Input: 
%   chnkr - chunker object
%
% Output:
%   ier - error flag
%       ier = 1, bad adj info, different number of left and right ends
%       ier = 2, bad adj info, missed/doubled chunks found
%
% Examples:
%   ier = checkadjinfo(chnkr)
%
% see also SORT, SORTINFO

% author: Travis Askham (askhamwhat@gmail.com)

[~,~,info] = sortinfo(chnkr);
ier = info.ier;