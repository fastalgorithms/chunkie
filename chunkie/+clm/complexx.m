
function [r,d,d2] = complexx(t,varargin)
%complexx
% return position, first and second derivatives of parameterized
% complexified layered media interface on the x-axis from [-L,L].
% 
%
% Inputs:
% t - parameter values to evaluate these quantities
%
% Optional inputs:
% r - coordinates
% d - first derivatives w.r.t. t
% d2 - second derivatives w.r.t. t

if nargin > 1 && ~isempty(varargin{1})
    L = varargin{1};
end
if nargin > 2 && ~isempty(varargin{2})
    c = varargin{2};
end
if nargin > 3 && ~isempty(varargin{3})
    a = varargin{3};
end

xs = t-c*1i*(erfc((t+L)/a)-erfc((L-t)/a));
ys = zeros(size(t));

xp = 1+c*1i*2/sqrt(pi)/a*(exp(-(t+L).^2/a^2)+exp(-(L-t).^2/a^2));
yp = ys;

xpp = -c*1i*4/sqrt(pi)/a^3*((t+L).*exp(-(t+L).^2/a^2)+(t-L).*exp(-(L-t).^2/a^2));
ypp = ys;

r = [(xs(:)).'; (ys(:)).'];
d = [(xp(:)).'; (yp(:)).'];
d2 = [(xpp(:)).'; (ypp(:)).'];

end

