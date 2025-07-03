function obj = helm2d_shape_der(type, zk, vfun)
%KERNEL.helm2d_shape_der  Construct the shape derivative matrices 
%   for helmholtz layer potentials when the perturbation is described
%   with a vector field stored as a data variable.
%
%   KERNEL.HELM2D_SHAPE_DER('s', zk, [datainds]) or 
%   KERNEL.HELM2D_SHAPE_DER('single', zk, [datainds]) 
%   constructs the shape derivate of the Helmholtz single-layer kernel.
%
%   KERNEL.HELM2D_SHAPE_DER('sp', zk, [datainds]) or 
%   KERNEL.HELM2D_SHAPE_DER('sprime', zk, [datainds]) 
%   constructs the shape derivate of the normal derivative of the
%   Helmholtz single-layer kernel.
%
% See also CHNK.HELM2D.KERN.

% author: Jeremy Hoskins

if nargin < 1
    error('Missing Helmholtz shape derivative kernel type.');
elseif nargin < 2
    error('Missing Helmholtz shape derivative kernel wavenumber.');
elseif nargin < 3 
    error('Missing data indices for Helmholtz shape derivative kernel.')
end
obj = kernel();
obj.name = 'helmholtz_shape_der';
obj.opdims = [1 1];
obj.params.zk = zk;
obj.params.datainds = vfun;

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.eval = @(s,t) chnk.helm2d_shape_der.kern(s, t, 's', zk, vfun);
        obj.sing = 'log';

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.eval = @(s,t) chnk.helm2d_shape_der.kern(s, t, 'sp', zk, vfun);
        obj.sing = 'log';


end

obj.fmm = [];

end


