function [r,d,d2] = sinecurve(tt,cpars)
%%circulararc
% return position, first and second derivatives of a sine curve
% that passes through points (x,y)=(a,0) and (x,y)=(b,0) with 
% amplitude A and number of half periods n
%
% Inputs:
% t - paramter values (-theta0/2,theta0/2) to evaluate these quantities
%
% Outputs:
% r - coordinates
% d - first derivatives w.r.t. t
% d2 - second derivatives w.r.t. t

a = cpars.a; b = cpars.b;
A = cpars.A; n = cpars.n;

lambda0 = n*pi;

lambda = 1;

islocal = -1;
if isfield(cpars,'islocal')
  islocal = cpars.islocal; 
end

if islocal == -1

  xs = a+(b-a)/lambda0*tt;
  ys = A*sin(lambda*tt);

  xp = (b-a)/lambda0*ones(size(tt));
  yp = A*lambda*cos(lambda*tt);

  xpp = zeros(size(tt));
  ypp = -lambda^2*ys;

elseif islocal == 1 % use stable and accurate forms to evaluate functions
% when the argument is close to 0. Note: series expansions are not always
% the best option.
  t = lambda*tt;
  st = sin(t);

  xs = (b-a)/lambda0*tt;
  ys = A*st;

  xp = (b-a)/lambda0*ones(size(tt));
  yp = A*lambda*cos(t);

  xpp = zeros(size(tt));
  ypp =  -lambda^2*ys;
elseif islocal == 0 
  t = lambda*tt;
  st = sin(t);
  
  mon = (-1)^n;
  xs = (b-a)/lambda0*tt;
  ys = A*mon*st;
  
  xp = (b-a)/lambda0*ones(size(tt));
  yp = A*mon*lambda*cos(t);

  xpp = zeros(size(tt));
  ypp =  -mon*lambda^2*ys;
  
end    
    

r = [(xs(:)).'; (ys(:)).'];
d = [(xp(:)).'; (yp(:)).'];
d2 = [(xpp(:)).'; (ypp(:)).'];

end