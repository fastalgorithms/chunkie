function obj = eye(dvals)
%KERNEL.EYE   Construct an identity ("Dirac delta") kernel, (K*sigma)(x) =
%   D(x)*sigma(x).
%
%   KERNEL.EYE() constructs the 1x1 identity delta.
%
%   KERNEL.EYE(A) constructs a constant delta with block A, a numeric
%   scalar or m x n matrix.
%
%   KERNEL.EYE(H) constructs h(t)*delta(t-s), H returning the block in
%   (m*nt x n) or (m x n x nt) format.
%
%   KERNEL.EYE(ZEROS(M,N)) can be used to make an operator with no identity
%   term (i.e. a first-kind equation).
%
%   See also KERNEL, KERNEL/PLUS, KERNEL/TIMES.

if nargin < 1 || isempty(dvals)
    dvals = 1;
end

if isa(dvals, 'function_handle')
    D0 = probe_handle(dvals);
    m = size(D0,1); n = size(D0,2);
    dfun = @(t) coerce_diag(dvals(t), m, n, size(t.r(:,:),2));
elseif isnumeric(dvals)
    if isscalar(dvals)
        m = 1; n = 1;
        A = dvals;
    else
        assert(ismatrix(dvals), 'KERNEL.EYE: matrix argument must be 2D.');
        m = size(dvals,1); n = size(dvals,2);
        A = dvals;
    end
    dfun = @(t) repmat(A, size(t.r(:,:),2), 1);
else
    error('KERNEL.EYE: argument must be a scalar, matrix, or function handle.');
end

obj        = kernel();
obj.name   = 'eye';
obj.type   = 'delta';
obj.opdims = [m n];
obj.sing   = 'smooth';
obj.iszero = false;
obj.isnan  = false;
obj.diag   = dfun;

obj.eval         = @(s,t) zeros(m*size(t.r(:,:),2), n*size(s.r(:,:),2));
obj.shifted_eval = @(s,t,o) zeros(m*size(t.r(:,:),2), n*size(s.r(:,:),2));

obj.fmm = @fmm_;
    function varargout = fmm_(eps, s, t, sigma) %#ok<INUSD>
        if ( isstruct(t) )
            nt = size(t.r(:,:), 2);
        else
            nt = size(t, 2);
        end
        if nargout > 0, varargout{1} = zeros(m*nt, 1); end
        if nargout > 1, varargout{2} = zeros(2, m*nt); end
        if nargout > 2, varargout{3} = zeros(3, m*nt); end
        if nargout > 3
            error('CHUNKIE:kernel:eye', 'Too many output arguments for FMM.');
        end
    end

end

function S = coerce_diag(D, m, n, nt)

if ismatrix(D) && size(D,1) == m*nt && size(D,2) == n
    S = D;
elseif size(D,1) == m && size(D,2) == n && size(D,3) == nt
    S = reshape(permute(D, [1 3 2]), m*nt, n);
else
    error(['KERNEL.EYE: diag handle must return (%d*nt x %d) stacked ', ...
           '(or [%d x %d x nt]); got size [%s].'], ...
           m, n, m, n, num2str(size(D)));
end
end

function D0 = probe_handle(h)

D0 = [];
for d = [2 3]
    try
        p.r = randn(d,1); p.n = randn(d,1); p.d = randn(d,1); p.d2 = randn(d,1);
        D0 = h(p);
        if ~isempty(D0), return; end
    catch
    end
end
if isempty(D0)
    error('KERNEL.EYE: unable to probe opdims from diag function handle.');
end
end
