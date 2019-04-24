function [val, der, der2] = absconvgauss(x, a, b, h)
%ABSCONVGAUSS this routine computes the convolution%
%  ( a*abs(x)+b ) \star 1/(sqrt(2*pi)*h^2) exp(-x^2/(2*h^2))
%
% this effectively smoothes off the corner from the abs(x) function
%
%

x2=x/sqrt(2.0)/h;
verf = erf(x2);
val=a*x.*verf+sqrt(2.0/pi)*a*h*exp(-x.*x/2.0/h/h)+b;
%
if nargout > 1; der=a*verf; end
if nargout > 2; der2=a*sqrt(2.0/pi)/h*exp(-x.*x/2.0/h/h); end
%
end
