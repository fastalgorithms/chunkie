function [gval, gdz, gdr, gdrp, gdrpr, gdzz, gdrz, gdrpz] = gfunc(r, rp, dr, z, zp, dz)

%
% chnk.axissymlap2d.gfunc evaluates the zeroth order axisymmetric Laplace
% Green's funcion, defined by the expression:
%
%     gfunc = pi * rp * \int_0^{2\pi} 1/|x - x'| d\theta'
%
% The extra factor of rp out front makes subsequent interfacing with RCIP
% slightly easier
%
% 
    domega3 = 2*pi;
    n       = 3;
    
    t = (dz.^2+dr.^2)./(2.*r.*rp);
    [q0,q1,q0d] = chnk.axissymlap2d.qleg_half(t);
    
    chi = t+1;
    %q0d = (-chi.*q0+q1)./(2*(chi.^2-1));
    
    gval = domega3*(rp./r).^((n-2)/2).*q0;
    gdz  = domega3*(rp./r).^((n-2)/2).*q0d ...
          ./(rp.*r).*dz;
    
    rfac = -r*(n-2)/2.*q0+(-(1+t).*r+rp).*q0d;
    gdrp  = domega3*(rp./r).^((n-2)/2)./(rp.*r) ...
            .*rfac;
    
    rfac = -rp*(n-2)/2.*q0+(-(1+t).*rp+r).*q0d;
    gdr  = domega3*(rp./r).^((n-2)/2)./(rp.*r) ...
           .*rfac;

    % Amandin hessian stuff
    chi = 1+t;
    %q0dd = ((1/4)*q0 + 2*chi.*q0d)./(1-chi.^2);
    q0dd = ((1/4)*q0 + 2*chi.*q0d)./(-t.^2-2*t);

    rfac = 1./(rp.*r).*q0d + (dz./(rp.*r)).^2.*q0dd;
    gdzz = 2*pi*sqrt(rp./r).*rfac;
    
    rfac = -3./(2*r.^2.*rp).*q0d + (-chi./(r.^2.*rp) + 1./(r.*rp.^2)).*q0dd;
    gdrz = 2*pi*sqrt(rp./r).*dz.*rfac;

    rfac = -3./(2*rp.^2.*r).*q0d + (-chi./(rp.^2.*r) + 1./(rp.*r.^2)).*q0dd;
    gdrpz = 2*pi*sqrt(rp./r).*dz.*rfac;

    rfac = 1./(4*r.*rp).*q0 ...
        + (2*chi./(rp.*r) - 3./(2*r.^2) - 3./(2*rp.^2)).*q0d ...
        + (-chi./r+1./rp).*(-chi./rp + 1./r).*q0dd;
    gdrpr = 2*pi*sqrt(rp./r).*rfac;
end