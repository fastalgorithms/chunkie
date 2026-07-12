function du = arclengthder(cgrph,u)
%ARCLENGTHDER arc length derivative of function on chunkgraph
%
% Calculates the arclength derivative of u using Gauss-Legendre 
%   spectral differentiation matrix. The derivative on the ith chunk is 
%   given by du(:,i) where du is the output of this routine
%
% Syntax: du = arclengthder(cgrph,u)
%
% Input:
%   cgrph - chunkgraph object
%   u - function or density defined on cgrph
%
% Output:
%   du - arclength derivative of u
%
% Examples:
%   du = arclengthder(cgrph,u);
%

dmat = lege.dermat(cgrph.k);

if cgrph.dim == 2
    ds = squeeze(sqrt(cgrph.d(1,:,:).^2+cgrph.d(2,:,:).^2));
elseif cgrph.dim == 3 
    ds = squeeze(sqrt(cgrph.d(1,:,:).^2+cgrph.d(2,:,:).^2+cgrph.d(3,:,:).^2));
end

du = dmat*reshape(u,[cgrph.k cgrph.nch]);
du = du ./ ds;

end
