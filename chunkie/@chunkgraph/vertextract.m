function [inds,isgn] = vertextract(ivert,cgrph)
%VERTEXTRACT for a single vertex in a chunkgraph object, sort the edges
% connecting to that vertex in counterclockwise order. The output is a 
% vector of edge numbers in sorted order, and a vector of signs indicating
% whether the given edge is incoming (plus in the incidence matrix) denoted
% by 1, or outgoing (minus in the incidence matrix) denoted by -1.
%
% Syntax: [inds,isgn] = vertextract(ivert,cgrph);
%
% Input:
%   ivert    - the index of the vertex to be processed
%   cgrph - chunkgraph object
%
% Output:
%   inds - a vector of edge numbers indicating which edges in cgrph 
%          connect to the vertex ivert. The elements of inds are sorted 
%          so that the corresponding edges are in counterclockwise order
%          around the ivert vertex.
%   isgn - a vector containing +1 in the ith entry if the ith edge in 
%          inds has its right end at ivert, and a -1 if the ith edge in 
%          inds has its left end at ivert.
%  
%
%

% author: Jeremy Hoskins

irel = find(cgrph.v2emat(:,ivert) ~= 0);
if isempty(irel)
     % if vertex is not connected, do nothing
    inds = [];
    isgn = [];
    return
end

% extract the indices of the edges which terminate at ivert.
ieplus = find(cgrph.edgesendverts(2,irel) ==  ivert);
ieplus = irel(ieplus); ieplus = ieplus(:).';
% extract the indices of the edges which begin at ivert.
ieminus  = find(cgrph.edgesendverts(1,irel) == ivert);
ieminus = irel(ieminus); ieminus = ieminus(:).';


% for each incoming edge, get the tangent vector near the end (at the 
% last discretization node)
deplus = [];
for i=1:numel(ieplus)
deplus = [deplus,-cgrph.echnks(ieplus(i)).d(:,end,end)];
end

% for each outgoing edge, get the tangent vector near the beginning (at the 
% first discretization node)
deminus = [];
for i=1:numel(ieminus)
deminus = [deminus,cgrph.echnks(ieminus(i)).d(:,1,1)];
end

% compute the corresponding angles and sort them
ds = [deplus,deminus];
angs = atan2(ds(2,:),ds(1,:));
[angs,isrtededges] = sort(angs);

% compute inds and isgn using isrtedges
inds = [ieplus,ieminus];
isgn = [ones(size(ieplus)),-ones(size(ieminus))];
inds = inds(isrtededges);
isgn = isgn(isrtededges);
end
