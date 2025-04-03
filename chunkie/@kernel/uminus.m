function f = uminus(f)
% multiply a kernel class by negative 1
%

if(isa(f.shifted_eval, 'function_handle'))
    f.shifted_eval = @(varargin) -f.shifted_eval(varargin{:});
else
    f.shifted_eval = [];

end

f.eval = @(varargin) - f.eval(varargin{:});

if (isa(f.fmm,'function_handle'))
    f.fmm = @(varargin) -f.fmm(varargin{:});
else
    f.fmm = [];
end

end


