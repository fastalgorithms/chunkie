function quiver(obj,varargin)
%QUIVER quiver plot of chunkgraph normal vectors in 2 dimensions
% Uses standard quiver commands. All edges in the chunkgraph will be
% plotted in a different color.
%
% Syntax: quiver(cgrph,varargin)
%
% Input: 
%   cgrph - chunkgraph object
%   varargin - any of the standard quiver commands
%
% Output:
%   none 
%
% Examples:
%   quiver(cgrph,'r','LineWidth',2) % quiver plot of chunkgraph normals
%                                   % with thick red arrows
%
% see also PLOT, PLOT3

% author: Jeremy Hoskins

ifhold = ishold();

hold on

for i=1:numel(obj.echnks)
    quiver(obj.echnks(i),varargin{:});
end    
    
hold off

if ifhold
    hold on
end
