function D = diffmat(chnkr, p, type)
%DIFFMAT spectral differentiation matrix with respect to arc length along a
% chunker object
%
% Syntax: chnkr = DIFFMAT(chnkr, p, type)
%
% Input:
%   chnkr - the chunker object
%
% Optional inputs:
%   p    - differentiation order (1)
%   type - type of spectral differentiation ('square')
%          - 'square' maps from an n-point grid to an n-point grid on each
%            panel.
%          - 'rectangular' maps from an n-point grid to an (n-p)-point grid
%            on each panel.
%
% Output:
%   D - chnkr.npt x chnkr.npt sparse block-diagonal differentiation matrix
%
% Examples:
%   D = diffmat(chnkr);
%   D = diffmat(chnkr, 2);
%   D = diffmat(chnkr, 'rect');
%
% See also INTMAT.

if ( nargin < 2 )
    p = 1;
    type = 'square';
end

if ( nargin < 3 )
    if ( isstring(p) || ischar(p) )
        type = p;
        p = 1;
    else
        type = 'square';
    end
end

if ( ~(isnumeric(p) && isscalar(p) && p>=0 && floor(p) == p) )
    error('CHUNKER:DIFFMAT:order', ...
        'Differentiation order must be a nonnegative integer.');
end

n = chnkr.k;
Dleg = lege.dermat(n);
ds = reshape(vecnorm(chnkr.d), [n 1 chnkr.nch]);
D = Dleg ./ ds;

% Take p-th derivative
for k = 1:chnkr.nch
    D(:,:,k) = D(:,:,k)^p;
end

if ( strcmpi(type, 'rect') )
    [~, ~, CfV, ~] = lege.exps(n);
    [~, ~, ~, VfC] = lege.exps(n-p);
	P = VfC * CfV(1:n-p,:);
	D = pagemtimes(P, D);
end

% Output a sparse block-diagonal matrix
D = num2cell(D, [1 2]);
D = matlab.internal.math.blkdiag(D{:});

end
