function obj = reverse(obj)
%REVERSE reverses the orientation of the chunker object
%
% Syntax: chnkr = reverse(chnkr)
%
% Input:
%   chnkr - chunker object
%
% Output:
%   chnkr - modified chunker object
%
% Examples:
%   chnkr = reverse(chnkr)
%
% see also SORT

% author: Travis Askham (askhamwhat@gmail.com)

k = obj.k;
obj.r = obj.r(:,k:-1:1,:);
obj.d = -obj.d(:,k:-1:1,:);
obj.d2 = obj.d2(:,k:-1:1,:);
obj.adj = obj.adj([2 1],:);
obj.n = -obj.n(:,k:-1:1,:);
obj.wts = obj.wts(k:-1:1,:);


end