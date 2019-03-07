
function [Kpxy,nbr] = pxyfun_circ(kern,proxy,pnorm,x,xnorm,slf,nbr,l,ctr)
%PXYFUN_CIRC proxy surface utility for circle proxy surface
%
% kern should be a kernel function of the form 
% submat = kern(src, targ, srcn, targn, slf)
%
% src and targ are (2,_) arrays of points
% srcn and targn are corresponding normals (same format)
% slf is the indices of the subset of the sources to use
%
% proxy should be a reference set of proxy points
% pnorm should be corresponding outward normals
% 
% x are the points of the geometry
% xnorm are the corresponding normals
%
% slf is a set of indices passed by the routine (relevant sources)
% nbr is a set of candidate points inside the proxy surface
% l length of a box at this level of the tree
% ctr center of this box
%

  pxy = bsxfun(@plus,proxy*l,ctr(:));
  N = size(x,2);
  Kpxy = kern(x,pxy,xnorm,pnorm,slf)*1.0/N;
  dx = x(1,nbr) - ctr(1);
  dy = x(2,nbr) - ctr(2);
  dist = sqrt(dx.^2 + dy.^2);
  nbr = nbr(dist/l < 1.5);
end
