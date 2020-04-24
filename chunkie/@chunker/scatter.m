function scatter(obj,varargin)
%SCATTER scatter plot of the coordinates of the chunker object. 
% Must be a 2D chunker. Faster than plot --> no sorting required
% Uses standard scatter commands
%
% Syntax: scatter(chnkr,varargin)
%
% Input:
%   chnkr - chunker object
%   varargin - any of the scatter commands
%
% Output:
%   none
%
% Examples:
%   scatter(chnkr,'ro') % scatter plot of chunker points as red circles
%
% see also PLOT

% author: Travis Askham (askhamwhat@gmail.com)

assert(obj.dim == 2,'for scatter plot must be 2D chunker');
xs = obj.r(1,:,:); xs = xs(:);
ys = obj.r(2,:,:); ys = ys(:);

scatter(xs,ys,varargin{:});