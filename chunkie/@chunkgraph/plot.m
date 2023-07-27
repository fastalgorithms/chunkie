function plot(obj,varargin)
%Plot plots chunkgraph in 2 dimensions
% Uses standard plot commands. All edges in the chunkgraph are
% plotted in a different color.
%
% Syntax: plot(cgrph,varargin)
%
% Input: 
%   cgrph - chunkgraph object
%   varargin - any of the standard plot commands
%
% Output:
%   none 
%
% Examples:
%   plot(cgrph,'r','LineWidth',2) % plot of chunkgraph normals
%                                   % with thick red arrows


% author: Jeremy Hoskins

ifhold = ishold();

echnks =  obj.echnks;

hold on
for i=1:numel(echnks)
    plot(echnks(i),varargin{:});
end

hold off

if ifhold
    hold on
end