function submat = kern(srcinfo, targinfo, origin, type)
%CHNK.AXISSYMLAP2D.KERN axissymmetric Laplace layer potential kernels in 2D
% 
% Syntax: submat = chnk.axissymlap2d.kern(srcinfo,targingo,type)
%
% Let x be targets and y be sources for these formulas, with
% n_x and n_y the corresponding unit normals at those points
% (if defined). Note that the normal information is obtained
% by taking the perpendicular to the provided tangential deriviative
% info and normalizing  
%
% Here the first and second components correspond to the r and z
% coordinates respectively. 
%
% Kernels based on G(x,y) = \int_{0}^{\pi} 1/(d(t)) \, dt \, 
% where d(t) = \sqrt(r^2 + r'^2 - 2rr' \cos(t) + (z-z')^2) with
% x = (r,z), and y = (r',z')
%
% D(x,y) = \nabla_{n_y} G(x,y)
% S(x,y) = G(x,y)
% S'(x,y) = \nabla_{n_x} G(x,y)
% D'(x,y) = \nabla_{n_x} \nabla_{n_y} G(x,y)
%
% Input:
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
%
%
% Output:
%   submat - the evaluation of the selected kernel for the
%            provided sources and targets. the number of
%            rows equals the number of targets and the
%            number of columns equals the number of sources  
%
%

src = srcinfo.r; 
targ = targinfo.r;

[~, ns] = size(src);
[~, nt] = size(targ);

if strcmpi(type, 'd')
    srcnorm = srcinfo.n;
    [~, grad] = chnk.axissymlap2d.green(src, targ, origin);
    nx = repmat(srcnorm(1,:), nt, 1);
    ny = repmat(srcnorm(2,:), nt, 1);
    % Due to lack of translation invariance in r, no sign flip needed, 
    % as gradient is computed with repsect to r'
    submat = (grad(:,:,2).*nx + grad(:,:,3).*ny);
end

if strcmpi(type, 'sprime')
    targnorm = targinfo.n;
    [~, grad] = chnk.axissymlap2d.green(src, targ, origin);

    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);
    submat = (grad(:,:,1).*nx - grad(:,:,3).*ny);

end

if strcmpi(type, 's')
    submat = chnk.axissymlap2d.green(src, targ, origin);
    
end


