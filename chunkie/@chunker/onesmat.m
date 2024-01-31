function mat = onesmat(chnkr)
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

wts = chnkr.wts; 
wts = wts(:);
temp = ones(size(wts));
mat = bsxfun(@times,temp,wts.');

end
