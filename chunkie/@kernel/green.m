function obj = green(green, type, varargin)
%KERNEL.GREEN   Construct a kernel from a Green's function.
%   KERNEL.GREEN(GREEN, TYPE) constructs the kernel of the given TYPE from
%   the Green's function G(x,y) = g(x-y). GREEN is a function handle in the
%   format of CHNK.HELM2D.GREEN, [val, grad, hess, ...] = GREEN(src, targ),
%   with derivatives with respect to the target. The types are 's', 'd',
%   'sp', 'dp', 'c', 'cp', 'sg', 'dg', 'cg', 'st', 'hilb', 's2trans',
%   'd2trans', 'c2trans', 'trans', 'trans_rep', 'trans_rep_prime',
%   'trans_rep_grad' and the plate kernels of CHNK.FLEX2D.KERN, see
%   CHNK.GREEN.KERN.
%
%   KERNEL.GREEN(GREEN, TYPE, COEFS) sets the coefficients of the combined
%   and transmission types, or the Poisson ratio of the plate kernels.
%
%   KERNEL.GREEN(GREEN, TYPE, ..., 'fmm', FMM) evaluates the kernel with
%   FMM(EPS, SRCINFO, PG, TARG, PGT), an FMM for GREEN with the interface
%   of HFMM2D without the wavenumber, see CHNK.GREEN.FMM.
%
%   KERNEL.GREEN(GREEN, TYPE, ..., 'dim', 3) is for points in 3D, and
%   'sing' sets the singularity type, which is otherwise found from the
%   kernel. 'check', false skips the check of the derivatives of GREEN.
%
%   KERNEL.GREEN(GREEN, FUN, 'nderiv', ND, 'opdims', [M N]) constructs a
%   custom kernel from the derivatives of GREEN, see CHNK.GREEN.KERN.
%
%   Sums and multiples of green kernels of the same built-in type are green
%   kernels of the sum of the Green's functions, with the sum of their fmms.
%
%   Example:
%      green = @(s,t) chnk.helm2d.green(zk, s, t);
%      fmm = @(eps,s,pg,t,pgt) hfmm2d(eps, zk, s, pg, t, pgt);
%      K = kernel.green(green, 'c', [1 1i], 'fmm', fmm);
%
% See also KERNEL.RADPP, CHNK.GREEN.KERN, CHNK.GREEN.FMM.

% kernel arithmetic: kernel.green('summable', KS) and
% kernel.green('sum', KS, CS) = sum_k CS(k)*KS{k}
if ( ischar(green) )
    if ( strcmp(green, 'summable') )
        obj = summable(type);
    else
        obj = gsum(type, varargin{1});
    end
    return
end

if ( nargin < 1 || ~isa(green, 'function_handle') )
    error('KERNEL.GREEN: the Green''s function must be a function handle.');
end
if ( nargin < 2 || isempty(type) )
    type = 's';
end

o = parseopts(varargin);
dim = o.dim;

if ( isa(type, 'function_handle') )
    if ( isempty(o.nderiv) || ~any(o.nderiv == 0:4) )
        error('KERNEL.GREEN: custom kernels need nderiv, 0 to 4.');
    end
    info = struct('type', 'custom', 'opdims', o.opdims, ...
        'nderiv', o.nderiv, 'param', '');
    ktype = struct('type', 'custom', 'fun', type, 'nderiv', o.nderiv, ...
        'opdims', o.opdims);
else
    info = chnk.green.typeinfo(type, dim);
    ktype = info;
end

param = typeparam(info, o);

if ( o.check )
    checkgreen(green, info.nderiv, dim);
end

obj = kernel();
obj.name = 'green';
obj.type = info.type;
obj.opdims = info.opdims;
obj.params.green = green;
obj.params.type = type;
obj.params.args = varargin;
obj.params.dim = dim;
obj.params.param = param;
obj.params.fmm = o.fmm;
obj.eval = @(s,t) chnk.green.kern(green, dim, s, t, ktype, param);

obj.fmm = [];
if ( ~isempty(o.fmm) && dim == 2 && chnk.green.fmm(info.type) )
    fmmfun = o.fmm;
    obj.fmm = @(eps,s,t,sigma) chnk.green.fmm(fmmfun, eps, s, t, ...
        info.type, sigma, param);
end

if ( isempty(o.sing) )
    obj.sing = singtype(obj.eval, info.opdims, dim);
else
    obj.sing = o.sing;
end

end


function tf = summable(Ks)
% green kernels of one built-in type, parameter and dimension
tf = all(cellfun(@(K) strcmp(K.name, 'green'), Ks));
if ( ~tf ), return; end
p = Ks{1}.params;
tf = all(cellfun(@(K) ~strcmp(K.type, 'custom') && strcmp(K.type, Ks{1}.type) && ...
    isequal(K.params.param, p.param) && K.params.dim == p.dim, Ks));
end


function obj = gsum(Ks, cs)
% sum_k cs(k)*Ks{k}, the kernel of the sum of the Green's functions
p = Ks{1}.params;
args = p.args;
i = find(cellfun(@(a) (ischar(a) || isstring(a)) && ...
    any(strcmpi(a, {'sing', 'fmm', 'check'})), args));
args([i, i+1]) = [];
sing = Ks{1}.sing;
for k = 2:numel(Ks), sing = kernel.singpromotion(sing, Ks{k}.sing); end
greens = cellfun(@(K) K.params.green, Ks, 'UniformOutput', false);
fmm = {};
if ( all(cellfun(@(K) ~isempty(K.params.fmm), Ks)) )
    fmms = cellfun(@(K) K.params.fmm, Ks, 'UniformOutput', false);
    fmm = {'fmm', chnk.green.sum(cs, fmms{:})};
end
obj = kernel.green(chnk.green.sum(cs, greens{:}), p.type, args{:}, fmm{:}, ...
    'sing', sing, 'check', false);
% otherwise the sum of the fmms of the kernels, if they have them
if ( isempty(fmm) && all(cellfun(@(K) isa(K.fmm, 'function_handle'), Ks)) )
    fmms = cellfun(@(K) K.fmm, Ks, 'UniformOutput', false);
    obj.fmm = chnk.green.sum(cs, fmms{:});
end
end


function o = parseopts(args)
% the positional parameter (coefs or nu) and the name-value options

o = struct('dim', 2, 'coefs', [], 'nu', [], 'sing', [], 'fmm', [], ...
    'check', true, 'nderiv', [], 'opdims', [1 1], 'par', []);

if ( ~isempty(args) && ~ischar(args{1}) && ~isstring(args{1}) )
    o.par = args{1};
    args(1) = [];
end
if ( mod(numel(args), 2) ~= 0 || ...
        ~all(cellfun(@(a) ischar(a) || isstring(a), args(1:2:end))) )
    error('KERNEL.GREEN: one parameter (coefs or nu) after the type, then name-value pairs.');
end
for i = 1:2:numel(args)
    name = lower(char(args{i}));
    if ( ~isfield(o, name) || strcmp(name, 'par') )
        error('KERNEL.GREEN: unknown option ''%s''.', args{i});
    end
    o.(name) = args{i+1};
end
if ( ~any(o.dim == [2 3]) )
    error('KERNEL.GREEN: dim must be 2 or 3.');
end

end


function param = typeparam(info, o)
% the parameter of the type: none, coefficients or nu

par = o.par;
switch info.param
    case ''
        if ( ~isempty(par) || ~isempty(o.coefs) || ~isempty(o.nu) )
            error('KERNEL.GREEN: type ''%s'' takes no coefs or nu.', info.type);
        end
        param = [];
    case 'nu'
        param = o.nu;
        if ( ~isempty(par) ), param = par; end
        if ( isempty(param) )
            error('KERNEL.GREEN: type ''%s'' needs nu.', info.type);
        end
    case 'coefs'
        param = o.coefs;
        if ( ~isempty(par) ), param = par; end
        if ( isempty(param) )
            warning('Missing coefficients. Defaulting to [1 1].');
            param = [1 1];
        end
        param = param(:).';
    case 'coefs22'
        param = o.coefs;
        if ( ~isempty(par) ), param = par; end
        if ( isempty(param) )
            warning('Missing coefficients. Defaulting to [1,1;1,1].');
            param = ones(2);
        end
        if ( numel(param) == 2 ), param = repmat(param(:).', 2, 1); end
end

end


function checkgreen(green, nderiv, dim)
% errors unless green returns nderiv+1 outputs of the right sizes and
% each derivative matches finite differences of the one below

err0 = 'KERNEL.GREEN: ';
names = {'val', 'grad', 'hess', 'third', 'fourth'};

nt = 3;  ns = 2;
s = randn(dim, ns);  t = randn(dim, nt) + [2; zeros(dim-1, 1)];
out = cell(1, nderiv+1);
try
    [out{:}] = green(s, t);
catch ME
    error([err0 'the type needs [%s] from green: %s'], ...
        strjoin(names(1:nderiv+1), ', '), ME.message);
end
if ( ~isequal(size(out{1}), [nt ns]) )
    error([err0 'val of green should be %s, not %s.'], ...
        mat2str([nt ns]), mat2str(size(out{1})));
end
for k = 1:nderiv
    try
        c = chnk.green.unpack(out{k+1}, dim, k);
    catch ME
        error([err0 '%s of green: %s'], names{k+1}, ME.message);
    end
    if ( size(c, 1) ~= nt*ns )
        error([err0 '%s of green should have %d points, not %d.'], ...
            names{k+1}, nt*ns, size(c, 1));
    end
end
if ( nderiv == 0 ), return; end

% each derivative against centred differences of the one below
y0 = [0.2; -0.1; 0.4];  y = y0(1:dim);
x0 = [0.5; 0.6; -0.3];  x = y + x0(1:dim);
h = 1e-5;
base = comps(green, y, x, nderiv, dim);
for k = 0:nderiv-1
    tup = tuples(dim, k);
    tup1 = tuples(dim, k+1);
    A = zeros(dim, size(tup, 1));
    F = zeros(dim, size(tup, 1));
    for j = 1:dim
        e = zeros(dim, 1);  e(j) = h;
        cp = comps(green, y, x+e, k, dim);
        cm = comps(green, y, x-e, k, dim);
        F(j,:) = (cp{k+1} - cm{k+1})/(2*h);
        for m = 1:size(tup, 1)
            A(j,m) = base{k+2}(all(tup1 == sort([tup(m,:) j]), 2));
        end
    end
    rel = norm(A - F, 'fro')/max(norm(A, 'fro'), realmin);
    if ( rel > 1e-6 )
        error([err0 '%s of green does not match finite differences ' ...
            'of %s (%.1e).'], names{k+2}, names{k+1}, rel);
    end
end

end


function tup = tuples(dim, k)
% nondecreasing index tuples of the derivative components of order k
tup = nchoosek(1:dim+k-1, k) - (0:k-1);
if ( k == 0 ), tup = zeros(1, 0); end
end


function C = comps(green, y, x, nderiv, dim)
% packed derivative components of green at one source and target
out = cell(1, nderiv+1);
[out{:}] = green(y, x);
C = cell(1, nderiv+1);
C{1} = out{1}(:).';
for k = 1:nderiv
    C{k+1} = chnk.green.unpack(out{k+1}, dim, k);
end
end


function sing = singtype(evalfun, opdims, dim)
% estimate the singularity type from the growth rho^-q of the kernel near
% a target on a parabola: q = 0 log, q = 1 pv, q >= 2 hs

rho = 10.^-(2:0.5:5);
nr = numel(rho);
x = [rho, -rho];
np = 2*nr;

ds = sqrt(1 + x.^2);
src = struct();
src.r = zeros(dim, np);   src.r(1,:) = x;   src.r(2,:) = x.^2/2;
src.d = zeros(dim, np);   src.d(1,:) = 1;   src.d(2,:) = x;
src.d2 = zeros(dim, np);  src.d2(2,:) = 1;
src.n = zeros(dim, np);   src.n(1,:) = x./ds;   src.n(2,:) = -1./ds;
src.data = 0.7*ones(2, np);

targ = struct();
targ.r = zeros(dim, 1);
targ.d = zeros(dim, 1);   targ.d(1) = 1;
targ.d2 = zeros(dim, 1);  targ.d2(2) = 1;
targ.n = zeros(dim, 1);   targ.n(2) = -1;
targ.data = 0.7*ones(2, 1);

try
    A = evalfun(src, targ);
    m = opdims(1);  n = opdims(2);
    mx = zeros(m*n, nr);
    for i = 1:m
        for j = 1:n
            v = abs(A(i,j:n:end));
            mx(i+m*(j-1),:) = max(v(1:nr), v(nr+1:end));
        end
    end
    if ( any(~isfinite(mx(:))) )
        error('the kernel is not finite on the test curve');
    end
    % entries below 1e-10 of the largest are zero up to rounding
    big = max(mx, [], 1);
    q = 0;
    for e = 1:m*n
        if ( all(mx(e,:) <= 1e-10*big) ), continue; end
        sl = diff(log(mx(e,:) + realmin))./diff(log(rho));
        q = max(q, max(0, ceil(-median(sl(end-2:end)) - 0.25)));
    end
catch ME
    warning('CHUNKIE:kernel:green:sing', ...
        'Could not determine the singularity of the kernel (%s). Using hs.', ...
        ME.message);
    q = 2;
end

sings = {'log', 'pv', 'hs'};
sing = sings{min(q, 2) + 1};

end
