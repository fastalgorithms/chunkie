function out = mtimes(f, g)
% * Matrix multiplication for kernel objects.
%
% Scalar c * K or K * c: scales eval, shifted_eval, fmm (same as times).
%
% Left multiply M * K: M is a (p x m) matrix or function handle M(t)
%   returning (p x m x nt). Output opdims = [p, K.opdims(2)].
%
% Right multiply K * N: N is a (q x p) matrix or function handle N(s)
%   returning (q x p x ns), q = K.opdims(2). Output opdims = [K.opdims(1), p].

% Determine which argument is the kernel and which is the multiplier.
if isa(f, 'kernel') && isa(g, 'kernel')
    error('KERNEL:mtimes:invalid', ...
        'Cannot * two kernel objects; use + to combine.');
end

if ~isa(f, 'kernel')
    [f, g] = deal(g, f);
    side = 'left';
elseif ~isa(g, 'kernel')
    side = 'right';
else
    error('KERNEL:mtimes:invalid', 'Unexpected argument types.');
end

K = f;
h = g;

% scalar: delegate to times
if isnumeric(h) && isscalar(h)
    out = times(K, h);
    return;
end

% constant matrix:
if isnumeric(h) && ~isscalar(h)
    A = h;
    if strcmp(side, 'left')
        assert(size(A,2) == K.opdims(1), ...
            'KERNEL:mtimes: left matrix must have %d columns', K.opdims(1));
        p = size(A, 1);
    else
        assert(size(A,1) == K.opdims(2), ...
            'KERNEL:mtimes: right matrix must have %d rows', K.opdims(2));
        p = size(A, 2);
    end
    h = A;
end

% only remaining option is a function handle
if ~isa(h, 'function_handle') && ~isnumeric(h)
    error('KERNEL:mtimes:invalid', ...
        'Argument must be a scalar, matrix, or function handle.');
end

if isa(h, 'function_handle')
    nargfunc = nargin(h);
    assert(nargfunc==1, ...
        'KERNEL:mtimes h must be a function of source or target, not both')
end

% probe h to determine output dimension p
%
% A function handle returning a (1 x 1 x n) array is treated as a
% pointwise scalar multiplier
ispointwise = false;
if ~exist('p', 'var')
    try
        probe.r  = randn(2,1);  probe.d  = randn(2,1);
        probe.d2 = randn(2,1);  probe.n  = randn(2,1);
        hval = h(probe);
        if size(hval,1) == 1 && size(hval,2) == 1
            ispointwise = true;
            if strcmp(side, 'left')
                p = K.opdims(1);
            else
                p = K.opdims(2);
            end
        elseif strcmp(side, 'left')
            p = size(hval, 1);
        else
            p = size(hval, 2);
        end
    catch
        error('KERNEL:mtimes:probe', ...
            'Could not probe function handle to determine output dimension.');
    end
end

if ~ispointwise && exist('hval', 'var')
    if strcmp(side, 'left')
        assert(size(hval,2) == K.opdims(1), ...
            'KERNEL:mtimes: left function handle must return matrices with %d columns', K.opdims(1));
    else
        assert(size(hval,1) == K.opdims(2), ...
            'KERNEL:mtimes: right function handle must return matrices with %d rows', K.opdims(2));
    end
end

Keval         = K.eval;
Kshifted_eval = K.shifted_eval;
Kfmm          = K.fmm;
m             = K.opdims(1);
q             = K.opdims(2);

out      = K;
out.type = ['custom_', K.type];
out.name = ['custom ', K.name];

    function fval = evalh(pts)
        % h is either a constant (p x m) or (q x p) matrix
        if isnumeric(h)
            fval = h;
        else
            fval = h(pts);
        end
    end

    function pts = shift_pts(pts, o)
        % Translate the positions of pts by o
        pts.r = pts.r + o(:);
    end

    function out = apply_left(fval, X)
        if size(fval,1) == 1 && size(fval,2) == 1
            out = fval .* X;
        else
            out = pagemtimes(fval, X);
        end
    end

    function out = apply_right(X, fval)
        if size(fval,1) == 1 && size(fval,2) == 1
            out = X .* fval;
        else
            out = pagemtimes(X, fval);
        end
    end

if strcmp(side, 'left')
    out.opdims       = [p, q];
    out.eval         = @eval_left;
    out.fmm          = set_if_exist(Kfmm, @fmm_left);
    out.shifted_eval = set_if_exist(Kshifted_eval, @shifted_eval_left);
else
    out.opdims       = [m, p];
    out.eval         = @eval_right;
    out.fmm          = set_if_exist(Kfmm, @fmm_right);
    out.shifted_eval = set_if_exist(Kshifted_eval, @shifted_eval_right);
end

% left-multiply: h(t) * K(s,t)

    function fval4 = reshape_fval_left(fval, nt)
        % fval is either p x m x nt (per-target) or p x m (constant);
        % reshape to p x m x nt x 1, broadcasting the constant case.
        if size(fval, 3) == nt
            fval4 = reshape(fval, size(fval,1), size(fval,2), nt, 1);
        else
            fval4 = reshape(fval, size(fval,1), size(fval,2), 1, 1);
        end
    end

    function vals = eval_left(s, t)
        nt   = size(t.r, 2);
        fval = evalh(t);
        Kmat = Keval(s, t);
        K4   = reshape(Kmat, m, 1, nt, []);
        out4 = apply_left(reshape_fval_left(fval, nt), K4);
        vals = reshape(out4, p*nt, []);
    end

    function vals = shifted_eval_left(s, t, o)
        nt   = size(t.r, 2);
        fval = evalh(shift_pts(t, o));
        Kmat = Kshifted_eval(s, t, o);
        K4   = reshape(Kmat, m, 1, nt, []);
        out4 = apply_left(reshape_fval_left(fval, nt), K4);
        vals = reshape(out4, p*nt, []);
    end

    function out = fmm_left(eps, s, t, sigma)
        nt    = size(t.r, 2);
        fval  = evalh(t);
        inner = Kfmm(eps, s, t, sigma);
        out   = reshape(apply_left(fval, reshape(inner, m, 1, nt)), p*nt, 1);
    end

% right-multiply: K(s,t) * h(s)

    function vals = eval_right(s, t)
        ns   = size(s.r, 2);
        nt   = size(t.r, 2);
        fval = evalh(s);
        Kmat = Keval(s, t);
        K3   = reshape(Kmat, m*nt, q, ns);
        vals = reshape(apply_right(K3, fval), m*nt, p*ns);
    end

    function vals = shifted_eval_right(s, t, o)
        ns   = size(s.r, 2);
        nt   = size(t.r, 2);
        fval = evalh(shift_pts(s, o));
        Kmat = Kshifted_eval(s, t, o);
        K3   = reshape(Kmat, m*nt, q, ns);
        vals = reshape(apply_right(K3, fval), m*nt, p*ns);
    end

    function out = fmm_right(eps, s, t, sigma)
        ns     = size(s.r, 2);
        fval   = evalh(s);
        sig_in = reshape(apply_left(fval, reshape(sigma, p, 1, ns)), q, ns);
        out    = Kfmm(eps, s, t, sig_in);
    end

end

function out = set_if_exist(cond, val)
if isa(cond, 'function_handle')
    out = val;
else
    out = [];
end
end
