function du = arclengthder(chnkr,u)
%ARCLENGTHDER arc length derivative of function on chunker
%
% Calculates the arclength derivative of u using Gauss-Legendre 
%   spectral differentiation matrix. The derivative on the ith chunk is 
%   given by du(:,i) where du is the output of this routine
%
% Syntax: du = arclengthder(chnkr,u)
%
% Input:
%   chnkr - chunker object
%   u - function or density defined on chnkr
%
% Output:
%   du - arclength derivative of u
%
% Examples:
%   du = arclengthder(chnkr,u);
%

dmat = lege.dermat(chnkr.k);
ds = squeeze(sqrt(chnkr.d(1,:,:).^2+chnkr.d(2,:,:).^2));

du = dmat*reshape(u,[chnkr.k chnkr.nch]);
du = du ./ ds;

end
