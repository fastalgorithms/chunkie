function f = conj(f)
% take a complex conjugate of a kernel class
%

if(isa(f.shifted_eval, 'function_handle'))
    f.shifted_eval = @(varargin) conj(f.shifted_eval(varargin{:}));
else
    f.shifted_eval = [];

end

f.eval = @(varargin) conj(f.eval(varargin{:}));

if (isa(f.fmm,'function_handle'))
    f.fmm = @(varargin) conj(f.fmm(varargin{:}));
else
    f.fmm = [];
end

end


