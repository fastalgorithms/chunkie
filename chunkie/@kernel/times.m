function f = times(f,g)
% .* Pointwise (elementwise) multiplication for kernel class
%
% K.*A and A.*K (equivalent) where A is numeric and K is a kernel with
% opdims [m q]. The product is taken elementwise on each source-target
% pair. A may be
%   - a scalar,
%   - a column vector (m x 1),
%   - a row vector (1 x q),
%   - an m x q matrix,
%
% Singleton dimensions of either A or K.opdims are expanded, just like the
% usual .* for matrices.
%
% The fmm adjusted when A is a scalar or column vector or a row vector.
% The fmm is set to [] if A is a matrix.

if (~isa(f,'kernel'))
    f = times(g,f);
    return
end

if isa(g,'kernel')
    error('KERNEL:times:invalid', ...
       'Cannot .* two kernel objects');
end

if ~isnumeric(g)
    error('KERNEL:times:invalid', ...
       'F or G must be numeric and the other a kernel class object');
end

if ~ismatrix(g)
    error('KERNEL:times:invalid', 'numeric factor must be a 2D array');
end

if isscalar(g)
    f = times_scalar(f,g);
    return
end

A = g;
m = f.opdims(1); q = f.opdims(2);
[a1, a2] = size(A);
if ~((a1 == m || a1 == 1 || m == 1) && (a2 == q || a2 == 1 || q == 1))
    error('KERNEL:times:dims', ...
        'KERNEL:times: size of A (%d x %d) incompatible with kernel opdims (%d x %d)', ...
        a1, a2, m, q);
end
P = max(a1, m); Q = max(a2, q);

if f.isnan || any(isnan(A(:)))
    f = kernel.nans(P, Q);
    return
end
if f.iszero || all(A(:) == 0)
    f = kernel.zeros(P, Q);
    return
end

Keval  = f.eval;
Kseval = f.shifted_eval;
Kfmm   = f.fmm;
A4 = reshape(A, a1, 1, a2, 1);

    function vals = apply_mat(Kmat)
        % Kmat is (m*nt) x (q*n); apply A blockwise with expansion
        nt = size(Kmat,1)/m;
        n  = size(Kmat,2)/q;
        K4 = reshape(Kmat, m, nt, q, n);
        vals = reshape(A4 .* K4, P*nt, Q*n);
    end

    function vals = eval_(s, t)
        vals = apply_mat(Keval(s, t));
    end

    function vals = shifted_eval_(s, t, o)
        vals = apply_mat(Kseval(s, t, o));
    end

    function out = fmm_col(eps, s, t, sigma)
        % A is a column vector: post-multiply the fmm output
        u  = Kfmm(eps, s, t, sigma);
        nt = numel(u)/m;
        out = reshape(A(:) .* reshape(u, m, nt), P*nt, 1);
    end

    function out = fmm_row(eps, s, t, sigma)
        % A is a row vector: pre-multiply the density
        ns  = numel(sigma)/Q;
        sig = A(:) .* reshape(sigma, Q, ns);
        if q == 1 && Q > 1
            % K has a single input channel: [a1 K, a2 K, ...] sigma
            sig = sum(sig, 1);
        end
        out = Kfmm(eps, s, t, reshape(sig, q*ns, 1));
    end

if isa(Keval, 'function_handle')
    f.eval = @eval_;
else
    f.eval = [];
end
if isa(Kseval, 'function_handle')
    f.shifted_eval = @shifted_eval_;
else
    f.shifted_eval = [];
end
if isa(Kfmm, 'function_handle') && a2 == 1
    f.fmm = @fmm_col;
elseif isa(Kfmm, 'function_handle') && a1 == 1
    f.fmm = @fmm_row;
else
    f.fmm = [];
end

s0 = f.splitinfo;
if isstruct(s0) && isfield(s0,'functions') && ~isempty(s0.functions) ...
        && isfield(s0,'type') && isfield(s0,'action')
    f0 = s0.functions;
    f.splitinfo = s0;
    f.splitinfo.functions = @(src,targ) cellfun(@apply_mat, ...
        f0(src,targ), 'UniformOutput', false);
else
    f.splitinfo = [];
end

f.opdims = [P, Q];
f.type = ['custom_', f.type];
f.name = ['custom ', f.name];

end

function f = times_scalar(f, g)
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

f.splitinfo = kernel.scale_splitinfo(f.splitinfo, g);

if or(f.isnan,isnan(g))
    f = kernel.nans(f.opdims(1),f.opdims(2));
end
if ~f.isnan && g==0 || f.iszero && ~isnan(g)
    f = kernel.zeros(f.opdims(1),f.opdims(2));
end
end
