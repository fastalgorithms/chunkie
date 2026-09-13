function f = mldivide(g, f)
% \ Left division for kernel objects
%
%   c \ K  where c is a numeric scalar        -> (1/c) * K
%   D \ K  where D is an m x m numeric matrix -> inv(D) * K
%   D \ K  where D is a function handle @(t)  -> inv(D(t)) * K
%
% m = K.opdims(1).

assert(isa(f,'kernel'), 'KERNEL:mldivide: denominator must be a kernel.');

if isnumeric(g) && isscalar(g)
    f = rdivide(f, g);
    return
end

m = f.opdims(1);

if isnumeric(g)
    assert(size(g,1) == size(g,2) && size(g,1) == m, ...
        'KERNEL:mldivide: matrix divisor must be %d x %d.', m, m);
    f = inv(g) * f;
    return
end

if isa(g,'function_handle')
    f = (@(t) pageinv_(g(t))) * f;
    return
end

error('KERNEL:mldivide: divisor must be a scalar, matrix, or function handle.');

end

function Ai = pageinv_(A)

if exist('pageinv','builtin') || exist('pageinv','file')
    Ai = pageinv(A);
    return
end
[m1, m2, nt] = size(A);
assert(m1 == m2, 'pageinv_: pages must be square.');
Ai = zeros(m1, m2, nt);
for k = 1:nt
    Ai(:,:,k) = inv(A(:,:,k));
end
end
