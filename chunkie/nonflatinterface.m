function [r,d,d2] = nonflatinterface(t,a,b,c,d)
%NONFLATINTERFACE curve parameterization of a perturbed interface given as 
% the graph (x(t),y(t)) = (t,d e^(-at^2/2) sin(bt+c)) 
% 
% Syntax: 
%    [r,d,d2] = nonflatinterface(t,a,b,c,d)
%
% Input: 
%    t - array of t parameter values where curve coordinates 
%    are desired
%    a,b,c,d - curve parameters as in formula above 
%
% Output:
%   r,d,d2 - 2 x numel(t) arrays of position, first and second 
%   derivatives at desired t parameter values
%

% author: Solomon Quinn 

   xs = t;
   xp = ones(size(t));
   xpp = zeros(size(t));
   
   
    ys = d*exp(-a*t.^2/2).*sin(b*t+c);
    yp = d*exp(-a*t.^2/2).*(b*cos(b*t+c) - a*t.*sin(b*t+c));
    ypp = d*exp(-a*t.^2/2).*(-2*a*b*t.*cos(b*t+c) + (a*a*t.^2-a-b^2).*sin(b*t+c));
    
    r = [(xs(:)).'; (ys(:)).'];
    d = [(xp(:)).'; (yp(:)).'];
    d2 = [(xpp(:)).'; (ypp(:)).'];
end
