function [f,fd,fdd] = complex_interface(t, width, slope, t0, t1)
%CHNK.CURVES.COMPLEX_INTERFACE evaluate the position, and first and second
% derivatives of a smooth interface curve defined via
%
%   f(t) = [t + i*(phi(t,a,b,t0) - phi(t,-a,b,t1)); 0]
%
% where a = width, b = slope/(2*width), and
% phi(t,u,v,z) = u*(t-z)*erfc(u*(t-z))*v - exp(-u^2*(t-z)^2)/sqrt(pi)*v
% is a smoothed piecewise linear function, that is below machine precision
% when width*(t-z) > 6.
%
% Syntax: [f,fd,fdd] = chnk.curves.complex_interface(t, width, slope, t0, t1)
%
% Input:
%   t     - (array) parameter values
%   width - controls sharpness of the transition (smaller = sharper)
%   slope - asymptotic slope of the imaginary part at large t
%   t0    - center of the first transition
%   t1    - center of the second transition

    a = width;
    b = slope / (2*width);
    phi   = @(t,u,v,z) u*(t-z).*erfc(u*(t-z))*v - exp(-u^2*(t-z).^2)/sqrt(pi)*v;
    phid  = @(t,u,v,z) u*erfc(u*(t-z))*v;
    phidd = @(t,u,v,z) -u*u*exp(-u^2*(t-z).^2)*2*v/sqrt(pi);
    f = zeros([2,size(t)]);
    fd = zeros([2,size(t)]);
    fdd = zeros([2,size(t)]);

    f(1,:) = t + 1i*(phi(t,a,b,t0) - phi(t,-a,b,t1));
    fd(1,:)= 1 + 1i*(phid(t,a,b,t0) - phid(t,-a,b,t1));
    fdd(1,:) = 1i*(phidd(t,a,b,t0) - phidd(t,-a,b,t1));

    f(2,:) = 0;
    fd(2,:) = 0;
    fdd(2,:) = 0;

end
