classdef kernel
%KERNEL class which describes an integral kernel, usually
% related to the solution of a partial differential equation (PDE)
% The kernel class should include a function, called kern,
% which takes as input two arguments (src,targ) in the ptinfo
% format. ptinfo is a struct with entries
%     src.r, src.d, src.d2, src.n
% corresponding to the position, first and second derivative,
% and normal to a point (the latter 3 are only defined for points
% on a parameterized curve).
%
% The class can also be used to store other useful methods
% and data for working with the kernel in a standardized
% format. E.g.:
%
% - a method called fastlayereval for rapidly evaluating
%  the integral over a chunker using only the smooth rule,
%  at a requested precision, e.g. by calling a compiled FMM
%  library
% - a method called kernsplits and a struct called splitinfo,
%  which provide the quantities needed to evaluate integrals
%  using the Helsing-Ojala technique  
%  
% author: Travis Askham (askhamwhat@gmail.com)

    properties(Access=private)
    end
    
    properties(Dependent,Access=public)
    end
    
    properties(Access=public)
      eval
      splitinfo
      fmm
      bdrysing
    end
    
    properties(Dependent,SetAccess=private)
    end
    
    methods
      function obj = kernel(varargin)
	if (nargin < 1)
	  obj = emptykernel();
	  return
	end
	v1 = varargin{1};
	if (isa(v1,'string') || isa(v1,'char'))
	  if (strcmpi(v1,'laplace') || strcmpi(v1,'l') || strcmpi(v1,'lap'))
	    obj = lap2dkernel(varargin{2:end});
	  elseif (strcmpi(v1,'zero'))
	    obj = zerokernel(varargin{2:end});
	  elseif (strcmpi(v1,'helmholtz') || strcmpi(v1,'helm') || strcmpi(v1,'h'))
	    obj = helm2dkernel(varargin{2:end});
	  elseif (strcmpi(v1,'stokes') || strcmpi(v1,'stok') || strcmpi(v1,'s'))
	    obj = stok2dkernel(varargin{2:end});
	  else
	    warning('kernel option (%s) not found, returning the zero kernel',v1);
	    obj = emptykernel();
	  end
	  return
	end
	if (isa(v1,'function_handle'))
	  obj = fun2kernel(varargin);
	  return
	end
	warning(['invalid first input in kernel constructor: ' ...
		   'first input should be a string or function ' ...
		   'handle. returning the zero kernel']);
	obj = emptykernel();
      end

      obj = lap2dkernel(type,opts);
      obj = emptykernel();
	
      
    end
end


function obj = emptykernel()
		     % return a default kernel
  obj.eval = @(s,t) [];
  obj.bdrysing = 'unknown';
end
