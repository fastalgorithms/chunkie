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

elseif islocal == 1 % series expansion around t=0
  if max(abs(tt))>0.1
    xs = (b-a)/lambda0*tt;
    ys = A*sin(lambda*tt);

    xp = (b-a)/lambda0*ones(size(tt));
    yp = A*lambda*cos(lambda*tt);

    xpp = zeros(size(tt));
    ypp = -lambda^2*ys;
  else
    t = lambda*tt;

    t2 = t.*t;
    t3 = t2.*t;
    t4 = t3.*t;
    t5 = t4.*t;
    t6 = t5.*t;
    t7 = t6.*t;
    t8 = t7.*t;
    t9 = t8.*t;
    t10 = t9.*t;  

    ct = -t2/2 + t4/24 - t6/720 + t8/40320 - t10/3628800;
    st =  t - t3/6 + t5/120 - t7/5040 + t9/362880;

    xs = (b-a)/lambda0*tt;
    ys = A*st;

    xp = (b-a)/lambda0*ones(size(tt));
    yp = A*lambda*(1+ct);

    xpp = zeros(size(tt));
    ypp =  -lambda^2*ys;
  end
elseif islocal == 0 % series expansion around t=1
  
  t = lambda*tt;
  
  t2 = t.*t;
  t3 = t2.*t;
  t4 = t3.*t;
  t5 = t4.*t;
  t6 = t5.*t;
  t7 = t6.*t;
  t8 = t7.*t;
  t9 = t8.*t;
  t10 = t9.*t;  

  ct = -t2/2 + t4/24 - t6/720 + t8/40320 - t10/3628800;
  st =  t - t3/6 + t5/120 - t7/5040 + t9/362880;
  
  mon = (-1)^n;
  xs = (b-a)/lambda0*tt;
  ys = A*mon*st;
  
  xp = (b-a)/lambda0*ones(size(tt));
  yp = A*mon*lambda*(1+ct);

  xpp = zeros(size(tt));
  ypp =  -mon*lambda^2*ys;
  
end    
    

r = [(xs(:)).'; (ys(:)).'];
d = [(xp(:)).'; (yp(:)).'];
d2 = [(xpp(:)).'; (ypp(:)).'];

end