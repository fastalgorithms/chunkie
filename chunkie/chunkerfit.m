function chnkr = chunkerfit(xy, opts)
%CHUNKERFIT   Create a chunker by fitting a curve to a set of points.
%
% Syntax: chnkr = CHUNKERFIT(pts, opts)
%
% Input:
%   xy - (2,:) array of points used to fit a curve
%
% Optional input:
%   opts - options structure (defaults)
%       opts.method = string ('spline')
%            Method to use for curve fitting.
%       opts.ifclosed = boolean (true)
%            Flag for fitting a closed (true) or open (false) curve.
%       opts.splitatpoints = boolean (false)
%            If true, start the adaptive splitting process from the set of
%            splits implied by the given points. This ensures that the given
%            points will always coincide with panel endpoints in the resulting
%            chunker object. If false, the adaptive splitting process will
%            start from scratch.
%       opts.cparams = struct
%            Curve parameters structure to be passed to chunkerfunc after
%            fitting is done. See chunkerfunc for parameters.
%       opts.pref = chunkerpref object or struct
%            Chunker preferences to be passed to chunkerfunc. See chunkerpref
%            for parameters.
%
% Output:
%   chnkr - chunker object discretizing the fitted curve
%
% Examples:
%   % Closed curve
%   r = chnk.curves.bymode(sort(2*pi*rand(20,1)), [2 0.5 0.2 0.7]);
%   chnkr = chunkerfit(r); % Closed curve
%   opts = []; opts.ifclosed = false;
%   chnkr = chunkerfit(r(:,1:10), opts); % Open curve
%
% See also CHUNKERFUNC, CHUNKERPREF, CHUNKER.

if ( nargin == 0 )
    test();
    return
end

if ( size(xy, 1) ~= 2 )
    error('CHUNKIE:CHUNKERFIT:size', 'Points must be specified as a 2xN matrix.');
end

x = xy(1,:).';
y = xy(2,:).';

% Set default options
defaults = [];
defaults.method = 'spline';
defaults.ifclosed  = true;
defaults.splitatpoints = false;
defaults.cparams = [];
defaults.pref = [];

if ( nargin < 2 )
    opts = defaults;
else
    for field = fieldnames(defaults).'
        if ( ~isfield(opts, field{1}) )
            opts.(field{1}) = defaults.(field{1});
        end
    end
end

if ( ~strcmpi(opts.method, 'spline') )
    error('CHUNKIE:CHUNKERFIT:method', 'Unsupported method ''%s''.', opts.method);
end

if ( opts.ifclosed )
    % Wrap the nodes around
    x  = [x; x(1)];
    y  = [y; y(1)];
end

% Choose the parametric grid to be polygonal arc length
t = [0; cumsum(sqrt(diff(x).^2 + diff(y).^2))];

ppx = myspline(t, x, opts.ifclosed);
ppy = myspline(t, y, opts.ifclosed);
ppx_dx  = ppdiff(ppx, 1);
ppy_dy  = ppdiff(ppy, 1);
ppx_dxx = ppdiff(ppx, 2);
ppy_dyy = ppdiff(ppy, 2);

% Pack polynomial coefficients into a tensor to vectorize ppval
cfs = zeros([6 size(ppx.coefs)]);
cfs(1,:,:) = ppx.coefs;
cfs(2,:,:) = ppy.coefs;
cfs(3,:,:) = ppx_dx.coefs;
cfs(4,:,:) = ppy_dy.coefs;
cfs(5,:,:) = ppx_dxx.coefs;
cfs(6,:,:) = ppy_dyy.coefs;
breaks = ppx.breaks.';

    function [r, d, d2] = splinefunc(tt)
        vals = myppval(breaks, cfs, tt);
        xs   = vals(1,:);
        ys   = vals(2,:);
        dxs  = vals(3,:);
        dys  = vals(4,:);
        d2xs = vals(5,:);
        d2ys = vals(6,:);
        r  = [xs(:).'  ; ys(:).'  ];
        d  = [dxs(:).' ; dys(:).' ];
        d2 = [d2xs(:).'; d2ys(:).'];
    end

cparams = opts.cparams;
cparams.ifclosed = opts.ifclosed;
cparams.ta = t(1);
cparams.tb = t(end);
if ( opts.splitatpoints )
    cparams.tsplits = t;
end
chnkr = chunkerfunc(@splinefunc, cparams, opts.pref);

end

function test()

rng(0)
r = chnk.curves.bymode(sort(2*pi*rand(20,1)), [2 0.5 0.2 0.7]);
chnkr = chunkerfit(r);

plot(chnkr, 'bo-')
hold on
plot(r(1,:), r(2,:), 'r.', markersize=30)
hold off
shg

end

function y = myppval(breaks, cfs, x)
    [~, idx] = histc(x, [-inf; breaks(2:end-1); inf]);
    p = reshape((x-breaks(idx)).^(3:-1:0), [1 size(idx,1) size(cfs,3)]);
    y = sum(p .* cfs(:,idx,:), 3);
end

function pp = ppdiff(pp, n)
    deg = pp.order-1;
    D = diag(deg:-1:1, 1);
    for k = 1:n
        pp.coefs = pp.coefs * D^n;
    end
end

function pp = myspline(x, y, ifclosed)
    n = length(y);
    dx = diff(x);
    dy = diff(y);
    dydx = dy ./ dx;
    b = zeros(n, 1);
    A = spdiags([  [dx(2:n-1) ; 0 ; 0],           ...
                 2*[0 ; dx(2:n-1)+dx(1:n-2) ; 0], ...
                   [0 ; 0 ; dx(1:n-2)] ],         ...
                [-1 0 1], n, n);

    if ( ifclosed )
        % Periodic conditions: match first and second derivatives at the first
        % and last points
        A(1,[1 n]) = [1 -1];
        A(n,1:2) = dx(n-1)*[2 1];
        A(n,n-1:n) =  A(n,n-1:n) + dx(1)*[1 2];
        b(2:n-1,:) = 3*(dx(2:n-1).*dydx(1:n-2) + dx(1:n-2).*dydx(2:n-1));
        b(n,:) = 3*(dx(n-1)*dydx(1) + dx(1)*dydx(n-1));
    else
        % Not-a-knot conditions: match third derivatives at the first and last
        % interior breaks
        A(1,1:2)   = [dx(2), x(3)-x(1)];
        A(n,n-1:n) = [x(n)-x(n-2), dx(n-2)];
        b(1,:) = ((dx(1)+2*(x(3)-x(1)))*dx(2)*dydx(1) + dx(1)^2*dydx(2)) / (x(3)-x(1));
        b(2:n-1,:) = 3*(dx(2:n-1).*dydx(1:n-2) + dx(1:n-2).*dydx(2:n-1));
        b(n,:) = (dx(n-1)^2*dydx(n-2) + ((2*(x(n)-x(n-2))+dx(n-1))*dx(n-2))*dydx(n-1)) / (x(n)-x(n-2));
    end

    % Solve linear system for the coefficients
    c = A \ b;
    c4 = (c(1:n-1)+c(2:n)-2*dydx(1:n-1)) ./ dx;
    c3 = (dydx(1:n-1)-c(1:n-1))./dx - c4;
    pp = mkpp(x, [c4./dx c3 c(1:n-1) y(1:n-1)]);
end
