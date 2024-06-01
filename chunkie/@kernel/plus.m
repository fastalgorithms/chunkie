function f = plus(f,g)
% + Pointwise addition for kernel class
%
% Currently only supported for adding two kernel class objects
% or adding a scalar to a kernel class object.

if (isa(g,'kernel') && isa(f,'kernel'))
  assert(f.opdims(1) == g.opdims(1) && f.opdims(1) == g.opdims(2), ...
      'kernel dimensions must agree to add');
  f.name = ['custom ',f.name,' ',g.name];
  f.eval = @(varargin) g.eval(varargin{:}) + f.eval(varargin{:});
  if (isa(g.fmm,'function_handle') && isa(f.fmm,'function_handle'))
    f.fmm = @(varargin) g.fmm(varargin{:}) + f.fmm(varargin{:});
  else
    f.fmm = [];
  end
  if or(f.isnan,g.isnan)
      f = kernel.nans(f.opdims(1),f.opdims(2));
  end
  if and(f.iszero,g.iszero)
      f.iszero = true;
  else
      f.iszero = false;
  end
else
    error('KERNEL:plus:invalid', ...
       'F and G must be either floats or kernel class objects');
end
end

