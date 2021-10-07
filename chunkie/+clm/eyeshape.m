function [r,d,d2] = eyeshape(t,icurve,cpars)
%%funcurve
% return position, first and second derivatives of a circular arc
% that passes through points (x,y)=(a,0) and (x,y)=(b,0) with opening
% angle theta0.
%
% Inputs:
% t - paramter values (-theta0/2,theta0/2) to evaluate these quantities
%
% Outputs:
% r - coordinates
% d - first derivatives w.r.t. t
% d2 - second derivatives w.r.t. t

a = cpars(1);b=cpars(2);theta0=cpars(3);

cx = (a+b)/2;
r0 = (b-a)/2/sin(theta0/2);

if icurve==3 % upper half
    cy = -(b-a)/2/tan(theta0/2);
    theta = pi/2+t;
elseif icurve==4 % bottom half
    cy = (b-a)/2/tan(theta0/2);
    theta = 3*pi/2+t;
end


xs = r0*cos(theta);
ys = r0*sin(theta);

xp = -ys;
yp = xs;

xpp = -xs;
ypp = -ys;

xs = cx + xs;
ys = cy + ys;


r = [(xs(:)).'; (ys(:)).'];
d = [(xp(:)).'; (yp(:)).'];
d2 = [(xpp(:)).'; (ypp(:)).'];

end