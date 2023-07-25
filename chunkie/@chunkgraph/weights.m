function wts = weights(chnkgrph)
%WEIGHTS integration weights suitable for smooth functions defined on the 
% chunkgraph object
%
% This is merely the standard Legendre weights scaled to the chunks
%
% Syntax: wts = weights(chnkgrph)
%
% Input:
%   chnkgrph - chunkgraph object
%
% Output:
%   wts - smooth integration weights
%
% Examples:
%   wts = weights(chnkgraph)
%

% author: Jeremy Hoskins


% merge the edge chunks to a single chunk object
  chnkr = merge(chnkgrph.echnks);
  
  k = chnkr.k;
  nch = chnkr.nch;
  % get the appropriate order Gauss-Legendre weights
  [~,w] = lege.exps(k);
  % get the norm of the tangent vector
  wts = reshape(sqrt(sum((chnkr.d).^2,1)),k,nch);
  wts = wts.*bsxfun(@times,w(:),(chnkr.h(:)).');

end
