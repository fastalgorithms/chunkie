classdef kernel
%KERNEL class which describes an integral kernel, usually
% related to the solution of a partial differential equation (PDE).
%
%   K = KERNEL(NAME, TYPE) constructs a kernel of the specified name and
%   type. The currently supported kernels names and types are:
%
%      NAME                           TYPE
%      ----                           ----
%      'laplace'    ('lap', 'l')      's', 'd', 'sp', 'c'
%      'helmholtz'  ('helm', 'h')     's', 'd', 'sp', 'c'
%      'elasticity' ('elast', 'e')    's', 'strac', 'd', 'dalt'
%      'stokes'     ('stok', 's')     'svel', 'spres', 'strac',
%                                     'dvel', 'dpres', 'dtrac'
%
%   The types may also be written in longer form, e.g. 'single', 'double',
%   'sprime', 'combined', 'svelocity', 'spressure', 'straction',
%   'dvelocity', 'dpressure', 'dtraction'.
%
%   Some kernels may require extra parameters to be specified, e.g. the
%   Helmholtz wavenumber, combined field parameter, or material parameters.
%   Such parameters should be passed as trailing arguments. For instance,
%   the Helmholtz single-layer kernel with wavenumber ZK can be constructed
%   as KERNEL('helmholtz', 's', ZK).
%
%   The kernel class includes an evaluator, K.eval(s,t), which evaluates
%   the kernel K between sources s and targets t specified in the ptinfo
%   format, i.e. a struct with entries s.r, s.d, s.d2, s.n corresponding to
%   the positions, first and second derivatives, and normals (the latter
%   three being only defined for points on a parameterized curve).
%
%   The kernel class can also be used to store other useful methods and
%   data for working with the kernel in a standardized format. For example,
%   for a kernel K:
%
%      - K.fmm: A function handle which calls the FMM for the corresponding
%        kernel. K.fmm(eps, s, t, sigma) evaluates the kernel with density
%        sigma from sources s to targets t with accuracy eps.
%
%      - K.splitinfo, which provides the quantities needed to evaluate
%        integrals using the Helsing-Ojala kernel-split quadrature
%        technique.

    properties

        name       % Name of the kernel
        type       % Type of the kernel
        params     % Structure of kernel parameters
        eval       % Function handle for kernel evaluation
        fmm        % Function handle for kernel FMM
        sing       % Singularity type
        splitinfo  % Kernel-split information

    end

    methods

      function obj = kernel(kern, varargin)

          if ( nargin < 1 )
              return
          end

          if ( isa(kern, 'string') || isa(kern, 'char') )
              switch lower(kern)
                  case {'laplace', 'lap', 'l'}
                      obj = kernel.lap2d(varargin{:});
                  case {'helmholtz', 'helm', 'h'}
                      obj = kernel.helm2d(varargin{:});
                  case {'stokes', 'stok', 's'}
                      obj = kernel.stok2d(varargin{:});
                  case {'elasticity', 'elast', 'e'}
                      obj = kernel.elast2d(varargin{:});
                  otherwise
                      error('Kernel ''%s'' not found.', kern);
              end
          elseif ( isa(kern, 'function_handle') )
              obj.eval = kern;
          else
              error('First input must be a string or function handle.');
          end

      end

    end

    methods ( Static )

        obj = lap2d(varargin);
        obj = helm2d(varargin);
        obj = stok2d(varargin);
        obj = elast2d(varargin);

    end

end
