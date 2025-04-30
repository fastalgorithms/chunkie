function f = mrdivide(f,g)
% * Matrix right division for kernel
%
% Currently only supports scalars
%
% returns c*F or F*C for the kernel F and scalar c
   f = rdivide(f,g);
end
