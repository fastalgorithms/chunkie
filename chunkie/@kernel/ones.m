function obj = ones(A)
%KERNEL.ONES   Construct a constant kernel.
%
%   K = KERNEL.ONES() constructs a 1 x 1 kernel with value 1, i.e.
%   K(x,y) = 1 for all targets x and sources y. Applied to a density this
%   gives the integral of the density: (K sigma)(x) = int sigma ds.
%
%   K = KERNEL.ONES(A) constructs an m x n constant block kernel with
%   K(x,y) = A, where A is an m x n numeric matrix.
%
%   A typical use is to remove the null space of a rank-deficient integral
%   equation, e.g. the interior Laplace Neumann problem
%
%       (1/2 I + S') sigma + W sigma = f,  W = kernel.ones()
%
%
%   See also KERNEL, KERNEL.ZEROS.

if ( nargin < 1 || isempty(A) )
    A = 1;
end
assert(isnumeric(A) && ismatrix(A), ...
    'CHUNKIE:kernel:ones', 'A must be a numeric 2D matrix.');

[m, n] = size(A);

    function out = eval_(s, t)
        % s.r / t.r may be 2 x k x nch (e.g. a chunker), so count points
        ns = size(s.r(:,:), 2);
        nt = size(t.r(:,:), 2);
        out = repmat(A, nt, ns);
    end

    function out = shifted_eval_(s, t, o) %#ok<INUSD>
        out = eval_(s, t);
    end

    function varargout = fmm_(eps, s, t, sigma) %#ok<INUSL>
        % sigma is the (quadrature-weighted) density, interleaved n x ns
        if ( isstruct(t) )
            nt = size(t.r(:,:), 2);
        else
            nt = size(t(:,:), 2);
        end
        blk = A * sum(reshape(sigma, n, []), 2);
        if ( nargout > 0 ), varargout{1} = repmat(blk, nt, 1); end
        if ( nargout > 1 ), varargout{2} = zeros(2, m*nt); end
        if ( nargout > 2 ), varargout{3} = zeros(3, m*nt); end
        if ( nargout > 3 )
            error('CHUNKIE:kernel:ones', 'Too many output arguments for FMM.');
        end
    end

obj = kernel();
obj.name         = 'ones';
obj.opdims       = [m n];
obj.sing         = 'smooth';
obj.eval         = @eval_;
obj.shifted_eval = @shifted_eval_;
obj.fmm          = @fmm_;
obj.iszero       = false;
obj.isnan        = false;

end
