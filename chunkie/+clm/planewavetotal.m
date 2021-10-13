function [u,gradu]=planewavetotal(k1,alpha,k2,targ,idomain,coef)

kstar = sqrt(k2^2-(k1*cos(alpha))^2);
c = coef(2)/coef(1)*k1*sin(alpha)/kstar;
R = (c-1)/(c+1);
%R = 2/(1+kstar*coef(2)/coef(1)/k1/sin(alpha))-1;
if idomain == 1
  x = targ(1,:,:); x = x(:);
  y = targ(2,:,:); y = y(:);
  t1 =   exp(1i*(k1*cos(alpha)*x-k1*sin(alpha)*y));
  t2 = R*exp(1i*(k1*cos(alpha)*x+k1*sin(alpha)*y));
  u = t1+t2;
  u = u(:);
  if nargout > 1
    gradu = zeros(length(x),2);
    gradu(:,1) =  1i*k1*cos(alpha)*u;
    gradu(:,2) = -1i*k1*sin(alpha)*(t1-t2);
  end
elseif idomain == 2
  x = targ(1,:,:); x = x(:);
  y = targ(2,:,:); y = y(:);
  u = (R+1)*exp(1i*(k1*cos(alpha)*x-kstar*y));
  u = u(:);
  if nargout > 1
    gradu = zeros(length(x),2);
    gradu(:,1) =  1i*k1*cos(alpha)*u;
    gradu(:,2) = -1i*kstar*u;
  end
end