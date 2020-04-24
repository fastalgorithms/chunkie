function mat = normonesmat(chnkr)
%NORMONESMAT for 2D curves. Forms the matrix for the operator 
%
%  W[\mu](x) = n(x) \int n(y)\cdot \mu(y) \, dl
%
% which project a vector density (\mu) on the chunker onto the normal
% vector along the chunker. This is a common operator in mechanics for
% removing nullspaces
%
% Syntax: mat = normonesmat(chnkr)
%
% Input:
%   chnkr - chunker object
%
% Output:
%   mat - the matrix discretization of the operator W above. mat is 
%       (2*chnkr.npt) x (2*chnkr.npt)
%
% Examples:
%   mat = normonesmat(chnkr)
%
% see also ONESMAT

% author: Travis Askham (askhamwhat@gmail.com)

wts = weights(chnkr);
rnorms = normals(chnkr);
wts = wts(:);
wts2 = repmat(wts.',2,1);
wts2 = wts2(:).*rnorms(:);

mat = bsxfun(@times,rnorms(:),wts2.');

end
