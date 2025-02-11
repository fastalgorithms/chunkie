function f = times(f,g)
% .* Pointwise multiplication for kernel class
%
% Currently only supported for scalars, if g is a scalar
% then f.*g or g.*f returns a kernel class object where each
% of the function handles - eval, and fmm are multiplied by 
% g

if (~isa(f,'kernel'))
    f = times(g,f);
    return
elseif (isscalar(g))
    if(isa(f.eval, 'function_handle'))        
        f.eval = @(varargin) g*f.eval(varargin{:});
    else
        f.eval = [];
    end
    
    if(isa(f.shifted_eval, 'function_handle'))        
        f.shifted_eval = @(varargin) g*f.shifted_eval(varargin{:});
    else
        f.shifted_eval = [];
    end
    
    if(isa(f.fmm, 'function_handle'))
        f.fmm = @(varargin) g*f.fmm(varargin{:});
    else
        f.fmm = [];
    end

    if or(f.isnan,isnan(g))
        f = kernel.nans(f.opdims(1),f.opdims(2));
    end
    if ~f.isnan && (g==0) || f.iszero && ~isnan(g)
        f = kernel.zeros(f.opdims(1),f.opdims(2));
    end
    
else
    error('KERNEL:times:invalid', ...
       'F or G must be scalar and the other a kernel class object');
end
end

