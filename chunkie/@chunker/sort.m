function [chnkr,info] = sort(chnkr)
%SORT sort the chunker object so that adjacent chunks have sequential 
% indices
% 
% Syntax: [chnkr,info] = sort(chnkr)
%
% Input:
%   chnkr - chunker object
%
% Output:
%   chnkr - the sorted chunker
%   info - adjacency info structure
%       info.ncomp - number of curve components detected
%       info.nchs - number of chunks on each component
%       info.ifclosed - whether/or not each detected component is closed
%       info.ier - error flag
%           ier = 1, bad adj info, different number of left and right ends
%           ier = 2, bad adj info, missed/doubled chunks found
%
% Examples:
%   chnkr = sort(chnkr)
%   [chnkr,info] = sort(chnkr)
%
% see also SORTINFO

% author: Travis Askham (askhamwhat@gmail.com)

[inds,adjs,info] = sortinfo(chnkr);

% reorder

chnkr.r = chnkr.r(:,:,inds);
chnkr.d = chnkr.d(:,:,inds);
chnkr.d2 = chnkr.d2(:,:,inds);
chnkr.h = chnkr.h(inds);
chnkr.n = chnkr.n(:,:,inds);

indinv = 1:length(inds);
indinv(inds) = 1:length(inds);

for i = 1:chnkr.nvert
    ivert = chnkr.vert{i};
    chnkr.vert{i} = indinv(ivert);
end

if chnkr.hasdata
    chnkr.data = chnkr.data(:,:,inds);
end

chnkr.adj = adjs;

    
end