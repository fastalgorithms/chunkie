function f = rdivide(f,g)
% ./ Pointwise division of kernel class
%
% Currently only supported for f a kernel and g a numeric scalar
%
% f./g returns a kernel class object where each of the function
% handles - eval, and fmm are divided by g

if (isnumeric(g) && isscalar(g))
    if(isa(f.eval, 'function_handle'))        
        f.eval = @(varargin) f.eval(varargin{:})/g;
    else
        f.eval = [];
    end
    
    if(isa(f.shifted_eval, 'function_handle'))        
        f.shifted_eval = @(varargin) f.shifted_eval(varargin{:})/g;
    else
        f.shifted_eval = [];
    end
    
    if(isa(f.fmm, 'function_handle'))
        f.fmm = @(varargin) f.fmm(varargin{:}).g;
    else
        f.fmm = [];
    end

    if or(f.isnan,isnan(g))
        f = kernel.nans(f.opdims(1),f.opdims(2));
    end
    if ~f.isnan && g==0 || f.iszero && ~isnan(g)
        f = kernel.zeros(f.opdims(1),f.opdims(2));
    end
    
else
    error('KERNEL:rdivide:invalid', ...
       'F must be a kernel class object and G a scalar');
end
end

