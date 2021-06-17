function nrm = normal2d(ptinfo)
%CHNK.NORMAL2D normal vector to curve for 2D curves
% 
% Syntax: nrm = chnk.normal2d(ptinfo)
%
% Input:
%   ptinfo - curve point info struct, with entries
%       ptinfo.r - positions (2,:) array
%       ptinfo.d - first derivative in underlying parameterization (2,:)
%       ptinfo.d2 - second derivative in underlying parameterization (2,:)
%
% Output:
%   nrm - (2,:) array containing corresponding normal information
%
% see also CHNK.CURVATURE2D

nrm = ptinfo.d;
dnrm = sqrt(sum(nrm.^2,1));
nrm = flipud(nrm);
nrm(2,:) = -nrm(2,:);
nrm = bsxfun(@rdivide,nrm,dnrm);