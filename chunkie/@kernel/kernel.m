classdef kernel
%KERNEL class which describes an integral kernel, usually
% related to the solution of a partial differential equation (PDE).
%
%   K = KERNEL(NAME, TYPE) constructs a kernel of the specified name and
%   type. The currently supported kernels names and types are:
%
%      NAME                                         TYPE
%      ----                                         ----
%      'laplace'    ('lap', 'l')                    's', 'd', 'sp', 'c'
%      'helmholtz'  ('helm', 'h')                   's', 'd', 'sp', 'dp', 'c'
%      'helmholtz1d'  ('helm1d', 'h1d')             's'
%      'helmholtz difference' ('helmdiff', 'hdiff') 's', 'd', 'sp', 'dp'
%      'elasticity' ('elast', 'e')                  's', 'strac', 'd', 'dalt'
%      'stokes'     ('stok', 's')                   'svel', 'spres', 'strac',
%                                                   'dvel', 'dpres', 'dtrac'
%      'zeros'       ('zero','z') 
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

        name           % Name of the kernel
        type           % Type of the kernel
        params         % Structure of kernel parameters
        eval           % Function handle for kernel evaluation
        fmm            % Function handle for kernel FMM
        sing           % Singularity type
        splitinfo      % Kernel-split information
        opdims = [0 0] % Dimension of the operator

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
                  case {'helmholtz1d', 'helm1d', 'h1d'}
                      obj = kernel.helm1d(varargin{:});
                  case {'helmholtz difference', 'helmdiff', 'hdiff'}
                      obj = kernel.helm2ddiff(varargin{:});
                  case {'stokes', 'stok', 's'}
                      obj = kernel.stok2d(varargin{:});
                  case {'elasticity', 'elast', 'e'}
                      obj = kernel.elast2d(varargin{:});
                  case {'zeros', 'zero', 'z'}
                      obj = kernel.zeros(varargin{:});
                  otherwise
                      error('Kernel ''%s'' not found.', kern);
              end
          elseif ( isa(kern, 'function_handle') )
              obj.eval = kern;
          elseif ( isa(kern, 'kernel') )
              if ( numel(kern) == 1 )
                  obj = kern;
              else
                  % The input is a matrix of kernels.
                  % Create a single kernel object by interleaving the
                  % outputs of each sub-kernel's eval() and fmm() routines.
                  % TODO: Check that opdims are consistent
                  obj = interleave(kern);
              end
          else
              error('First input must be a string or function handle.');
          end

      end

    end

    methods ( Static )

        obj = lap2d(varargin);
        obj = helm2d(varargin);
        obj = helm1d(varargin);
        obj = helm2ddiff(varargin);
        obj = stok2d(varargin);
        obj = elast2d(varargin);
        obj = zeros(varargin);

    end

end

function K = interleave(kerns)
%INTERLEAVE   Create a kernel from a matrix of kernels by interleaving.

[m, n] = size(kerns);

rowdims = cat(1, kerns(:,1).opdims); rowdims = rowdims(:,1).';
coldims = cat(1, kerns(1,:).opdims); coldims = coldims(:,2).';
rowstarts = [0 cumsum(rowdims)];
colstarts = [0 cumsum(coldims)];
opdims = [sum(rowdims) sum(coldims)];

    function out = eval_(s, t)

        [~, ns] = size(s.r);
        [~, nt] = size(t.r);

        % Compute interleaved indices
        ridx = cell(m, 1);
        cidx = cell(n, 1);
        for k = 1:m
            ridx{k} = (rowstarts(k)+1:rowstarts(k+1)).' + (0:nt-1)*opdims(1);
            ridx{k} = ridx{k}(:).';
        end
        for l = 1:n
            cidx{l} = ((colstarts(l)+1):colstarts(l+1)).' + (0:ns-1)*opdims(2);
            cidx{l} = cidx{l}(:).';
        end

        % Evaluate each sub-kernel and assign the resulting block to the
        % output matrix using interleaved indices
        out = zeros(opdims(1)*nt, opdims(2)*ns);
        for k = 1:m
            for l = 1:n
                out(ridx{k},cidx{l}) = kerns(k,l).eval(s,t);  
            end
        end

    end

    function varargout = fmm_(eps, s, t, sigma)

        [~, ns] = size(s.r);
        if isa(t,'struct')
            [~,nt] = size(t.r);
        else
            [~,nt] = size(t);
        end
        

        % Compute interleaved indices
        ridx = cell(m, 1);
        cidx = cell(n, 1);
        for k = 1:m
            ridx{k} = (rowstarts(k)+1:rowstarts(k+1)).' + (0:nt-1)*opdims(1);
            ridx{k} = ridx{k}(:).';
        end
        for l = 1:n
            cidx{l} = ((colstarts(l)+1):colstarts(l+1)).' + (0:ns-1)*opdims(2);
            cidx{l} = cidx{l}(:).';
        end

        % Call the FMM for each sub-kernel using interleaved indices to
        % slice the given density
        fmms = cell(m, n, nargout);
        for k = 1:m
            for l = 1:n
                [fmms{k,l,:}] = kerns(k,l).fmm(eps, s, t, sigma(cidx{l}));
            end
        end

        % Combine the outputs from each block using interleaved indices
        if ( nargout > 0 )
            pot = zeros(opdims(1)*nt, 1);
            for k = 1:m
                for l = 1:n
                    pot(ridx{k}) = pot(ridx{k}) + fmms{k,l,1};
                end
            end
            varargout{1} = pot;
        end

        if ( nargout > 1 )
            grad = zeros(2, opdims(1)*nt);
            for k = 1:m
                for l = 1:n
                    grad(:,ridx{k}) = grad(:,ridx{k}) + fmms{k,l,2};
                end
            end
            varargout{2} = grad;
        end

        if ( nargout > 2 )
            error('CHUNKIE:kernel:interleave', 'Too many output arguments for FMM.');
        end

    end

K = kernel();
K.opdims = opdims;

% The new singularity type is the worst of the singularity types of all
% sub-kernels
sings = {kerns.sing};
K.sing = 'smooth';
if ( any(strcmpi(sings, 'log')) ),  K.sing = 'log'; end
if ( any(strcmpi(sings, 'pv'))  ),  K.sing = 'pv';  end
if ( any(strcmpi(sings, 'hs'))  ),  K.sing = 'hs';  end

% The new kernel has eval() only if all sub-kernels have eval()
if ( all(cellfun('isclass', {kerns.eval}, 'function_handle')) )
    K.eval = @eval_;
end

% The new kernel has fmm() only if all sub-kernels have fmm()
if ( all(cellfun('isclass', {kerns.fmm}, 'function_handle')) )
    K.fmm  = @fmm_;
end

end
