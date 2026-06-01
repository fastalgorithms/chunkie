function [r, d, d2] = ellipse(t,varargin) 
%ELLIPSE return position, first and second derivatives of an ellipse 
% with the parameterization 
%
% x(t) = a*cos(t)
% y(t) = b*sin(t)
%
% Syntax: [r,d,d2] = ellipse(t,a,b) 
%
% Input:
%   t - array of points (in [0,2pi])
%
% Optional input:
%   a - xscaling
%   b - xscaling
%
% Output:
%   r - 2 x numel(t) array of positions, r(:,i) = [x(t(i)); y(t(i))]
%   d - 2 x numel(t) array of t derivative of r 
%   d2 - 2 x numel(t) array of second t derivative of r 
%
% Examples:
%   [r,d,d2] = ellipse(t); % circle parameterization
%   [r,d,d2] = ellipse(t,a,b); % stretch circle into ellipse
%
a = 1;
b = 1;
if nargin > 1 && ~isempty(varargin{1})
    a = varargin{1};
end
if nargin > 2 && ~isempty(varargin{2})
    b = varargin{2};
end

r =  [ a*cos(t(:).'); b*sin(t(:).')];
d =  [-a*sin(t(:).'); b*cos(t(:).')];
d2 = [-a*cos(t(:).');-b*sin(t(:).')];



end