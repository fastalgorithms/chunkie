function obj = radpp(green, rmax, type, varargin)
%KERNEL.RADPP   Construct a kernel from a piecewise polynomial fit of a
% radial Green's function.
%   KERNEL.RADPP(GREEN, RMAX, TYPE) constructs the kernel of the given TYPE
%   from a fit of the radial Green's function GREEN for distances up to
%   RMAX, or in RMAX = [RMIN RMAX] (default RMIN = 1e-12). GREEN, TYPE and
%   the options are as in KERNEL.GREEN. The kernel is NaN outside the fit.
%
%   KERNEL.RADPP(GREEN, RMAX, TYPE, ..., 'opts', OPTS) sets the options of
%   the fit, see CHNK.RADPP.FIT. The fit is stored in K.params.
%
%   Sums and multiples of radpp kernels of the same type are radpp kernels
%   whose fit is the sum of the fits, on the distances common to them.
%
%   Example:
%      green = @(s,t) chnk.helm2d.green(zk, s, t);
%      K = kernel.radpp(green, [1e-6 10], 'sp');
%
% See also KERNEL.GREEN, CHNK.RADPP.FIT, CHNK.RADPP.GREEN.

% kernel arithmetic: kernel.radpp('summable', KS) and
% kernel.radpp('sum', KS, CS) = sum_k CS(k)*KS{k}
if ( ischar(green) )
    if ( strcmp(green, 'summable') )
        obj = summable(rmax);
    else
        obj = ppsum(rmax, type);
    end
    return
end

if ( nargin < 2 )
    error('Missing Green''s function or range of distances.');
end
if ( nargin < 3 )
    type = 's';
end

% the fit options are ours, the rest are those of kernel.green
fitopts = struct();
args = varargin;
i = find(cellfun(@(a) (ischar(a) || isstring(a)) && strcmpi(a, 'opts'), args), 1);
if ( ~isempty(i) )
    fitopts = args{i+1};
    args(i:i+1) = [];
end
if ( ~isfield(fitopts, 'eps') ), fitopts.eps = 1e-13; end

% the exact kernel checks the inputs and gives the singularity type
Kex = kernel.green(green, type, args{:});
dim = Kex.params.dim;
if ( isa(type, 'function_handle') )
    nderiv = args{find(cellfun(@(a) ischar(a) && strcmpi(a, 'nderiv'), args), 1) + 1};
else
    info = chnk.green.typeinfo(type, dim);
    nderiv = info.nderiv;
end

% the rotated exact profiles must be green
exact = @(r) profiles(green, r, nderiv, dim);
Krot = kernel.green(@(s,t) chnk.radpp.green(exact, dim, s, t), type, ...
    args{:}, 'sing', Kex.sing, 'check', false);
checkradial(Kex, Krot, dim);

rmin = 1e-12;
if ( numel(rmax) == 2 )
    rmin = rmax(1);
    rmax = rmax(2);
end
fit = chnk.radpp.fit(exact, rmin, rmax, fitopts);

breaks = fit.breaks;
coefs = fit.coefs;
ppgreen = @(s,t) chnk.radpp.green( ...
    @(r) chnk.radpp.fiteval(r, breaks, coefs), dim, s, t);
obj = kernel.green(ppgreen, type, args{:}, 'sing', Kex.sing, 'check', false);

obj.name = 'radpp';
obj.params.green = green;
obj.params.ppgreen = ppgreen;
obj.params.type = type;
obj.params.args = args;
obj.params.fitopts = fitopts;
obj.params.rmin = rmin;
obj.params.rmax = rmax;
obj.params.breaks = breaks;
obj.params.coefs = coefs;
obj.params.nbin = fit.nbin;
obj.params.nlev = fit.nlev;
obj.params.err = fit.err;

end


function tf = summable(Ks)
% radpp kernels of one built-in type, parameter and fit order, with fits
% that overlap
tf = all(cellfun(@(K) strcmp(K.name, 'radpp'), Ks));
if ( ~tf ), return; end
p = Ks{1}.params;
tf = all(cellfun(@(K) ~strcmp(K.type, 'custom') && strcmp(K.type, Ks{1}.type) && ...
    isequal(K.params.param, p.param) && K.params.dim == p.dim && ...
    size(K.params.coefs, 2) == size(p.coefs, 2) && ...
    size(K.params.coefs, 3) == size(p.coefs, 3), Ks));
tf = tf && max(cellfun(@(K) K.params.rmin, Ks)) < ...
    min(cellfun(@(K) K.params.rmax, Ks));
end


function obj = ppsum(Ks, cs)
% sum_k cs(k)*Ks{k}, with the fits added on the union of their panels
p = Ks{1}.params;
rmin = max(cellfun(@(K) K.params.rmin, Ks));
rmax = min(cellfun(@(K) K.params.rmax, Ks));
breaks = unique(cell2mat(cellfun(@(K) K.params.breaks(:), Ks(:), ...
    'UniformOutput', false)));
breaks = unique([rmin; breaks(breaks > rmin & breaks < rmax); rmax]);
a = breaks(1:end-1);
[~, n, m] = size(p.coefs);
coefs = zeros(numel(a), n, m);
for k = 1:numel(Ks)
    % the panel polynomials, highest degree first in r - b, re-centred at
    % a: (r-b)^q = ((r-a) + (a-b))^q
    q = Ks{k}.params;
    ib = discretize((a + breaks(2:end))/2, q.breaks);
    d = q.coefs(ib, end:-1:1, :);
    del = a - q.breaks(ib);
    e = zeros(size(d));
    for j = 0:n-1
        for i = j:n-1
            e(:,j+1,:) = e(:,j+1,:) + nchoosek(i, j)*del.^(i-j).*d(:,i+1,:);
        end
    end
    coefs = coefs + cs(k)*e(:, end:-1:1, :);
end

sing = Ks{1}.sing;
for k = 2:numel(Ks), sing = kernel.singpromotion(sing, Ks{k}.sing); end
par = {};
if ( ~isempty(p.param) ), par = {p.param}; end
ppgreen = @(s,t) chnk.radpp.green( ...
    @(r) chnk.radpp.fiteval(r, breaks, coefs), p.dim, s, t);
% the fmm of the sum of the greens
fmm = {};
if ( all(cellfun(@(K) ~isempty(K.params.fmm), Ks)) )
    fmms = cellfun(@(K) K.params.fmm, Ks, 'UniformOutput', false);
    fmm = {'fmm', chnk.green.sum(cs, fmms{:})};
end
obj = kernel.green(ppgreen, Ks{1}.type, par{:}, 'dim', p.dim, fmm{:}, ...
    'sing', sing, 'check', false);
% otherwise the sum of the fmms of the kernels, if they have them
if ( isempty(fmm) && all(cellfun(@(K) isa(K.fmm, 'function_handle'), Ks)) )
    fmms = cellfun(@(K) K.fmm, Ks, 'UniformOutput', false);
    obj.fmm = chnk.green.sum(cs, fmms{:});
end

greens = cellfun(@(K) K.params.green, Ks, 'UniformOutput', false);
obj.name = 'radpp';
obj.params.green = chnk.green.sum(cs, greens{:});
obj.params.ppgreen = ppgreen;
obj.params.type = Ks{1}.type;
obj.params.args = [par, {'dim', p.dim}, fmm];
obj.params.fitopts = p.fitopts;
obj.params.rmin = rmin;
obj.params.rmax = rmax;
obj.params.breaks = breaks;
obj.params.coefs = coefs;
obj.params.nbin = numel(a);
obj.params.nlev = max(cellfun(@(K) K.params.nlev, Ks));
obj.params.err = sum(abs(cs(:)).*cellfun(@(K) K.params.err, Ks(:)));
end


function f = profiles(green, r, nderiv, dim)
% the on-axis profiles of green, (numel(r),nprof), with the source at the
% origin and the targets on the x-axis, see CHNK.RADPP.GREEN

% the components of the derivatives of each order with an even number of
% y's, in the lexicographic order of the packed derivatives
sel = {1, [1 1; 2 2], [1 1 1; 1 2 2], [1 1 1 1; 1 1 2 2; 2 2 2 2]};

r = r(:);
out = cell(1, nderiv+1);
[out{:}] = green(zeros(dim, 1), [r.'; zeros(dim-1, numel(r))]);

f = reshape(out{1}, [], 1);
for k = 1:nderiv
    c = chnk.green.unpack(out{k+1}, dim, k);
    tup = nchoosek(1:dim+k-1, k) - (0:k-1);
    for j = 1:size(sel{k}, 1)
        f = [f, c(:, all(tup == sel{k}(j,:), 2))]; %#ok<AGROW>
    end
end

end


function checkradial(Kex, Krot, dim)
% the kernel from the rotated on-axis profiles must be that of green
s = struct('r', randn(dim, 4), 'n', [], 'd', [], 'd2', [], 'data', []);
t = struct('r', randn(dim, 3) + [2; zeros(dim-1, 1)]);
s.d = randn(dim, 4);  s.d2 = randn(dim, 4);  s.data = randn(2, 4);
t.d = randn(dim, 3);  t.d2 = randn(dim, 3);  t.data = randn(2, 3);
s.n = s.d./vecnorm(s.d);  t.n = t.d./vecnorm(t.d);
A = Kex.eval(s, t);
B = Krot.eval(s, t);
if ( norm(A - B, 'fro') > 1e-8*norm(A, 'fro') )
    error('KERNEL.RADPP: green is not radial.');
end
end
