function obj = radpp(green, rmax, type, varargin)
%KERNEL.RADPP   Construct a kernel from a piecewise polynomial fit of a
% radial Green's function.
%   KERNEL.RADPP(GREEN, RMAX, TYPE) constructs the kernel of the given TYPE
%   for the radial Green's function G(|x-y|) given by the function handle
%   GREEN. The types, the format of GREEN and the options are as in
%   KERNEL.RADGREEN, which evaluates GREEN directly. KERNEL.RADPP does not
%   call GREEN when the kernel is evaluated.
%
%   Piecewise polynomial fit of the required derivatives:
%
%      [G | G_x | G_xx G_yy | G_xxx G_xyy | G_xxxx G_xxyy G_yyyy]
%
%   up to the order needed by the type, for distances in RMAX or
%   [RMIN RMAX]. The default RMIN is 1e-12. The kernel is NaN for
%   distances outside the fit.
%
%   KERNEL.RADPP(GREEN, RMAX, TYPE, ..., 'opts', OPTS) sets the options of
%   the fit, see CHNK.RADPP.FIT. The default tolerance OPTS.eps is 1e-13.
%
%   Custom kernels are as in KERNEL.RADGREEN.
%
%   The fit is stored in K.params.fit.
%
%   Example:
%      zk = 1.1;
%      green = @(s,t) chnk.helm2d.green(zk, s, t);
%      K = kernel.radpp(green, [1e-6 10], 'sp');
%
% See also KERNEL.RADGREEN, CHNK.RADPP.FIT, CHNK.RADPP.SUMGREENS.

if ( nargin < 2 )
    error('Missing Green''s function or range of distances.');
end
if ( nargin < 3 )
    type = 's';
end

obj = chnk.radpp.build(green, rmax, type, varargin, true);

end
