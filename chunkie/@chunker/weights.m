function wts = weights(chnkr)
%WEIGHTS compute integration weights suitable for smooth functions defined 
% on the chunker object.  This routine is intended for use by chunker object
% constructors. It is not intended to be called by users, because 
% such weights should already be stored in the chunker object. 
%
% This is the standard Legendre weights scaled to the chunks
%
% Syntax: wts = weights(chnkr)
%
% Input:
%   chnkr - chunker object
%
% Output:
%   wts - smooth integration weights
%
% Examples:
%   wts = weights(chnkr)
%

% author: Travis Askham (askhamwhat@gmail.com)

  k = chnkr.k;
  nch = chnkr.nch;
  w = chnkr.wstor;
  wts = reshape(sqrt(sum((chnkr.d).^2,1)),k,nch);
  wts = wts.*w(:);

end
