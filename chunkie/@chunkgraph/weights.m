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
  wts = chnkr.wts;

end
