function scatter(obj,varargin)
%SCATTER scatter plot of the coordinates of the chunkgraph object. 
% Must be a 2D chunkgraph. Faster than plot --> no sorting required
% Uses standard scatter commands
%
% Syntax: scatter(obj,varargin)
%
% Input:
%   obj - chunkgraph object
%   varargin - any of the scatter commands
%
% Output:
%   none
%
% Examples:
%   scatter(obj,'ro') % scatter plot of chunkgraph points as red circles
%
% see also PLOT

% author: Travis Askham (askhamwhat@gmail.com)

assert(obj.dim == 2,'for scatter plot must be 2D chunkgraph');
xs = obj.r(1,:,:); xs = xs(:);
ys = obj.r(2,:,:); ys = ys(:);

scatter(xs,ys,varargin{:});