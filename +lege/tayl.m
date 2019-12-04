function [polnew,dernew] = tayl(pol,der,x,h,n,k)
%
%        This subroutine evaluates the Taylor series for the
%        Legendre polynomial and its derivative at the point
%        x+h, starting with the value of the polynomial at the
%        point x, and the value of the derivative of that 
%        polynomial. It uses the obvious three-term recursion 
%        for the derivatives of Legendre polynomials.
% c
%                 Input parameters:
% c
%  pol - the value of the polynomial at the pount x
%  der - the derivative of the polynomial at the pount x
%  x - the point where the values pol, der are specified
%  h - the polynomial and its derivative will be evaluated at 
%        the point x+h
%  n - the order of the Legendre polynomial
%  k - the order of the Taylor series to be used
% c
%                 Output parameters:
% c
%  polnew - the value of P_n at the point x+h
%  dernew - the value of the derivative of P_n at the point x+h
% c
%        . . . evaluate the derivatives of P_n scaled by h^n/n!,
%              and polnew the taylor series for P_n and its 
%              derivative
% c

q0 = pol;
q1 = der.*h;
q2 = (2*x.*der-n*(n+1)*pol)/(1-x.^2);
q2 = q2.*h.^2/2.0;

polnew=q0+q1+q2;
dernew=q1./h+(q2*2)./h;

if(k <= 2); return; end
 
qi=q1;
qip1=q2;

for i = 1:k-2
    d=2*(x*(i+1).^2)./h.*qip1-(n*(n+1)-i*(i+1))*qi;
    d=(d/(i+1)/(i+2)).*h.^2./(1-x.^2);
    qip2=d;
    polnew=polnew+qip2;
    dernew=dernew+(d*(i+2))./h;
    qi=qip1;
    qip1=qip2;
end

end