function [r,d,d2] = linesegment(t,cpars)
%%circulararc
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

v0 = cpars.v0; v1 = cpars.v1;
x0 = v0(1); y0 = v0(2);
x1 = v1(1); y1 = v1(2);

islocal = -1;
if isfield(cpars,'islocal')
  islocal = cpars.islocal; 
end

if islocal == -1
  xs = x0 + (x1-x0)*t;
  ys = y0 + (y1-y0)*t;

  xp = (x1-x0)*ones(size(t));
  yp = (y1-y0)*ones(size(t));

  xpp = zeros(size(t));
  ypp = zeros(size(t));
else % local series expansion around t=0 or t=1, removing the constant term
  xs = (x1-x0)*t;
  ys = (y1-y0)*t;

  xp = (x1-x0)*ones(size(t));
  yp = (y1-y0)*ones(size(t));

  xpp = zeros(size(t));
  ypp = zeros(size(t));
end    
    
r = [(xs(:)).'; (ys(:)).'];
d = [(xp(:)).'; (yp(:)).'];
d2 = [(xpp(:)).'; (ypp(:)).'];

end