function obj = lap2d_shape_der(type, vfun)
%KERNEL.LAP2D_SHAPE_DER  Construct the shape derivative matrices 
%   for helmholtz layer potentials when the perturbation is described
%   with a vector field stored as a data variable.
%
%   KERNEL.LAP2D_SHAPE_DER('s', zk, [datainds]) or 
%   KERNEL.LAP2D_SHAPE_DER('single', zk, [datainds]) 
%   constructs the shape derivate of the Laplace single-layer kernel.
%
%   KERNEL.LAP2D_SHAPE_DER('sp', zk, [datainds]) or 
%   KERNEL.LAP2D_SHAPE_DER('sprime', zk, [datainds]) 
%   constructs the shape derivate of the normal derivative of the
%   Laplace single-layer kernel.
%
% See also CHNK.HELM2D.KERN.

% author: Jeremy Hoskins

if nargin < 1
    error('Missing Laplace shape derivative kernel type.');
elseif nargin < 2
    error('Missing data indices for Laplace shape derivative kernel.')
end


obj = kernel();
obj.name = 'laplace_shape_der';
obj.opdims = [1 1];
obj.params.datainds = vfun;

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.eval = @(s,t) chnk.lap2d_shape_der.kern(s, t, 's', vfun);
        obj.sing = 'log';

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.eval = @(s,t) chnk.lap2d_shape_der.kern(s, t, 'sp', vfun);
        obj.sing = 'log';

end

obj.fmm = [];

end


