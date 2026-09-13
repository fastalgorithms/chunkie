function f = rdivide(f,g)
% ./ Pointwise division of kernel class
%
% Currently only supported for f a kernel and g a numeric scalar

if (isnumeric(g) && isscalar(g))
    if g == 0
        f = times(f, 0);
    else
        f = times(f, 1/g);
    end
else
    error('KERNEL:rdivide:invalid', ...
       'F must be a kernel class object and G a scalar');
end
end
