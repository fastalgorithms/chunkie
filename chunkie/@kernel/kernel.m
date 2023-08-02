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
        opdims     % Dimension of the operator

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
          elseif ( isa(kern, 'kernel') )
              % TODO: Check that opdims are consistent
              kern = interleave(kern);
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

function K = interleave(kern)

[m, n] = size(kern);
K = kernel();

opdims = [0 0];

opdims_cat = reshape(cat(3,K.opdims),[1,2 size(K)]);
opdims_rows = squeeze(opdims_cat(1,1,1:m,1));
opdims_cols = squeeze(opdims_cat(1,2,1,1:n));
opdims_rows = opdims_rows(:);
opdims_cols = opdims_cols(:);

opdims(1) = sum(opdims_rows);
opdims(2) = sum(opdims_cols);

opdims_rows_csum = [0; cumsum(opdims_rows)]';
opdims_cols_csum = [0; cumsum(opdims_cols)]';


K.opdims = opdims;


    function K_interleaved = k_eval(s, t)
        [~, ns] = size(s.r);
        [~, nt] = size(t.r);
        
        irinds = cell(m,1);
        icinds = cell(n,1);
        
        for k=1:m
            irinds{k} = (opdims_rows_csum(k)+1):(opdims_rows_csum(k+1)) +  ...
               (0:(nt-1))*K.opdims(1);
        end
        for l=1:n
            icinds{l} = (opdims_cols_csum(l)+1):(opdims_rows_csum(l+1)) +  ...
               (0:(ns-1))*K.opdims(2);
        end    
             
        K_interleaved = zeros(nt*K.opdims(1), ns*K.opdims(2));
        for k = 1:m
            for l = 1:n
                K_interleaved(irinds{k},icinds{l}) = kern(k,l).eval(s,t);  
            end
        end
    end

    function varargout = k_fmm(eps, s, t, sigma)
        
        fmm_kl = cell(nargout, m, n);
        
        [~, ns] = size(s.r);
        [~, nt] = size(t.r);
        
        irinds = cell(m,1);
        icinds = cell(n,1);
        
        for k=1:m
            irinds{k} = (opdims_rows_csum(k)+1):(opdims_rows_csum(k+1)) +  ...
               (0:(nt-1))*K.opdims(1);
        end
        for l=1:n
            icinds{k} = (opdims_cols_csum(l)+1):(opdims_rows_csum(l+1)) +  ...
               (0:(ns-1))*K.opdims(2);
        end    
        
        for k = 1:m
            for l = 1:n
            [fmm_kl{:,k,l}] = kern(k,l).fmm(eps, s, t, sigma(icinds{l}));
            end
        end
        
        if(nargout >=1)
            pot = zeros(K.opdims(1)*nt,1);
            for k=1:m
                for l=1:n
                    pot(irinds{k}) = pot(irinds{k}) + fmm_kl{1,k,l};
                end
            end
            varargout{1} = pot;
        elseif(nargout >=2)
            grad = zeros(2,K.opdims(1)*nt);
            for k=1:m
                for l=1:n
                    grad(:,irinds{k}) = grad(:,irinds{k}) + fmm_kl{2,k,l};
                end
            end
            varargout{2} = grad;
        else
            error('In KERNEL.INTERLEAVE: Too many output arguments for fmm, aborting.\n');
        end   
    end

K.eval = @k_eval;
K.fmm = @k_fmm;

end
