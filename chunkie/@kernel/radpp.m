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
%   Example:
%      green = @(s,t) chnk.helm2d.green(zk, s, t);
%      K = kernel.radpp(green, [1e-6 10], 'sp');
%
% See also KERNEL.GREEN, CHNK.RADPP.FIT, CHNK.RADPP.GREEN.

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
obj.params.rmin = rmin;
obj.params.rmax = rmax;
obj.params.breaks = breaks;
obj.params.coefs = coefs;
obj.params.nbin = fit.nbin;
obj.params.nlev = fit.nlev;
obj.params.err = fit.err;

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
