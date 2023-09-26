
function submat= kern(zk,srcinfo,targinfo,type,varargin)
%CHNK.HELM1D.KERN standard Helmholtz layer potential kernels in 2D
% 
% Syntax: submat = chnk.heml1d.kern(zk,srcinfo,targingo,type,varargin)
%
% Let x be targets and y be sources for these formulas, with
% n_x and n_y the corresponding unit normals at those points
% (if defined). Note that the normal information is obtained
% by taking the perpendicular to the provided tangential deriviative
% info and normalizing  
%  
% Kernels based on G(x,y) =  e^{i zk  |t-t'|} where t parametrization
% coordinate
%
% D(x,y) = \nabla_{n_y} G(x,y)
% S(x,y) = G(x,y)
% S'(x,y) = \nabla_{n_x} G(x,y)
%
% Input:
%   zk - complex number, Helmholtz wave number
%   srcinfo - description of sources in ptinfo struct format, i.e.
%                ptinfo.r - positions (2,:) array
%                ptinfo.d - first derivative in underlying
%                     parameterization (2,:)
%                ptinfo.d2 - second derivative in underlying
%                     parameterization (2,:)
%   targinfo - description of targets in ptinfo struct format,
%                if info not relevant (d/d2) it doesn't need to
%                be provided. sprime requires tangent info in
%                targinfo.d
%   type - string, determines kernel type
%                type == 'd', double layer kernel D
%                type == 's', single layer kernel S
%                type == 'sprime', normal derivative of single
%                      layer S'
%                type == 'c', combined layer kernel D + i eta S
%   varargin{1} - eta in the combined layer formula, otherwise
%                does nothing
%
% Output:
%   submat - the evaluation of the selected kernel for the
%            provided sources and targets. the number of
%            rows equals the number of targets and the
%            number of columns equals the number of sources  
%
% see also CHNK.HELM2D.GREEN
  
src = srcinfo.r(1,:);
src = src(:).';
targ = targinfo.r(1,:);
targ = targ(:);

[~,ns] = size(src);
[~,nt] = size(targ);


if strcmpi(type,'s')
  submat = exp(1j*zk*abs(targ-src))/2/1j/zk;
end
end


