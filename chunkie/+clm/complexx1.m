
function [r,d,d2] = complexx1(t,icurve,cpars)
%complexx
% return position, first and second derivatives of parameterized
% complexified layered media interface on the x-axis from [-L,b]
% or [b,L].
%
% Inputs:
% t - parameter values to evaluate these quantities
%
% Optional inputs:
% r - coordinates
% d - first derivatives w.r.t. t
% d2 - second derivatives w.r.t. t

L = cpars(1);c1=cpars(2);c2=cpars(3);


if icurve==1 % left part
  
   xs = t-c1*1i*erfc((t+L)/c2);
   xp = 1+c1*1i*2/sqrt(pi)/c2*exp(-(t+L).^2/c2^2);
   xpp = -c1*1i*4/sqrt(pi)/c2^3*(t+L).*exp(-(t+L).^2/c2^2);
   
elseif icurve==2 % right part
  
   xs = t+c1*1i*erfc((L-t)/c2);
   xp = 1+c1*1i*2/sqrt(pi)/c2*exp(-(L-t).^2/c2^2);
   xpp = -c1*1i*4/sqrt(pi)/c2^3*(t-L).*exp(-(L-t).^2/c2^2);   
end

ys = zeros(size(t));
yp = ys;
ypp = ys;

r = [(xs(:)).'; (ys(:)).'];
d = [(xp(:)).'; (yp(:)).'];
d2 = [(xpp(:)).'; (ypp(:)).'];

end

