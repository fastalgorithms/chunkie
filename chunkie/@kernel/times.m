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
        feval_str = func2str(f.eval);  
        fstr_split = split(feval_str,')');
        fstr_split{2} = [num2str(double(g)) '*' fstr_split{2}];
        fstr_join = join(fstr_split,')');
        f.eval = str2func(fstr_join{1});
    else
        f.eval = [];
    end
    
    if(isa(f.fmm, 'function_handle'))
        ffmm_str = func2str(f.fmm);  
        fstr_split = split(ffmm_str,')');
        fstr_split{2} = [num2str(double(g)) '*' fstr_split{2}];
        fstr_join = join(fstr_split,')');
        f.fmm = str2func(fstr_join{1});
    else
        f.fmm = [];
    end
    
else
    error('KERNEL:times:invalid', ...
       'F or G must be scalar and the other a kernel class object');
end
end

