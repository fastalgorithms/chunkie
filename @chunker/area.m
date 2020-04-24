function a = area(chnkr)
%AREA compute area of a closed, 2D chunker curve
%
% Syntax: a = area(chnkr)
%
% Input:
%   chnkr - chunker object. Must have chnkr.dim==2, all(chnkr.adj) ~= 0,
%           i.e. must be a closed 2d curve (or collection of closed curves)
%
% Output:
%   a - area
%
% Examples:
%   a = area(chnkr)
%

% author: Travis Askham (askhamwhat@gmail.com)

assert(chnkr.dim == 2, ... 
    'area only well-defined for 2d chunkers');
assert(nnz(chnkr.adj == 0) == 0, ...
    'area only well-defined for closed 2d chunkers');
assert(~any(chnkr.vertdeg > 2), ...
    'area not well-defined for higher order vertices');

wts = weights(chnkr);
rnorm = normals(chnkr);
igrand = sum(rnorm.*(chnkr.r),1); igrand = igrand(:);
a = sum(wts(:).*igrand)/chnkr.dim;

end