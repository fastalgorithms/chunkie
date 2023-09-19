function ds = arclengthdens(chnkr)
%ARCLENGTHDENS arc length density on chunks
%
% warning: this takes the length of the ith chunk in parameter space to 
%   be 2*chnkr.h(i) as opposed to 2. Thus the smooth integration rule
%   on the ith chunk is 
%       ds(:,i).*w*chnkr.h(i)
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
