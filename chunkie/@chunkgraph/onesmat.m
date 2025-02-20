function mat = onesmat(cgrph)
%ONESMAT Forms the matrix for the operator 
%
%  W[\mu](x) = \int \mu(y) dl
%
% which projects a scalar density (\mu) on the chunker onto the constant
% vector along the chunker. This is a common operator for removing 
% nullspaces
%
% Syntax: mat = onesmat(chnkr)
%
% Input:
%   cgrph - chunkgraph object
%
% Output:
%   mat - the matrix discretization of the operator W above. mat is 
%       (cgrph.npt) x (cgrph.npt)
%
% Examples:
%   mat = onesmat(chnkr)
%
% see also NORMONESMAT

% author: Travis Askham (askhamwhat@gmail.com)

wts = cgrph.wts; 
wts = wts(:);
temp = ones(size(wts));
mat = bsxfun(@times,temp,wts.');

end
