function ds = arclengthdens(chnkr)
%ARCLENGTHDENS arc length density on chunks
%
%   The smooth integration rule
%   on the ith chunk is 
%       ds(:,i).*w
%   where w are the standard Legendre weights of appropriate order and 
%   ds is the output of this routine
%
% Syntax: ds = arclengthdens(chnkr)
%
% Input:
%   chnkr - chunker object
%
% Output:
%   ds - arclength density
%
% Examples:
%   ds = arclengthdens(chnkr);
%
% see also WEIGHTS

% author: Travis Askham (askhamwhat@gmail.com)

ds = squeeze(sqrt(sum((chnkr.d).^2,1)));

end
