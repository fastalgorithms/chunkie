function quiver(obj,varargin)
%QUIVER quiver plot of chunker normal vectors in 2 dimensions
% Uses standard quiver commands
%
% Syntax: quiver(chnkr,varargin)
%
% Input: 
%   chnkr - chunker object
%   varargin - any of the standard quiver commands
%
% Output:
%   none 
%
% Examples:
%   quiver(chnkr,'r','LineWidth',2) % quiver plot of chunker normals
%                                   % with thick red arrows
%
% see also PLOT, PLOT3

% author: Travis Askham (askhamwhat@gmail.com)

ifhold = ishold();

hold on

for i=1:numel(obj.echnks)
    quiver(obj.echnks(i),varargin{:});
end    
    
hold off

if ifhold
    hold on
end
