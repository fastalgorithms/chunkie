function [gval,gdz,gdr,gdrp] = gfunc(r, rp, dr, z, zp, dz)
%
% chnk.axissymlap2d.gfunc evaluates the zeroth order axisymmetric Laplace
% Green's funcion, defined by the expression:
%
%     gfunc = pi * rp * \int_0^{2\pi} 1/|x - x'| d\theta
%
% The extra factor of rp out front makes subsequent interfacing with RCIP
% slightly easier
%
% 
    domega3 = 2*pi;
    n       = 3;

    t = (dz.^2+dr.^2)./(2.*r.*rp);
    [q0,q1,q0d] = chnk.axissymlap2d.qleg_half(t);

    gval = domega3*(rp./r).^((n-2)/2).*q0;
    gdz  =-domega3*(rp./r).^((n-2)/2).*q0d ...
           ./(rp.*r).*dz;

    rfac = -r*(n-2)/2.*q0+(-(1+t).*r+rp).*q0d;
    gdrp  = -domega3*(rp./r).^((n-2)/2)./(rp.*r) ...
          .*rfac;

    rfac = -rp*(n-2)/2.*q0+(-(1+t).*rp+r).*q0d;
    gdr  = -domega3*(rp./r).^((n-2)/2)./(rp.*r) ...
          .*rfac;


    % % check the scaling here correctly
    % m = 1000;
    % dint = 0;
    % h = 2*pi/m;
    % for i = 1:m
    %     s = (i-1)*h;
    %     x = r;
    %     y = 0;
    %     xp = rp*cos(s);
    %     yp = rp*sin(s);
    %     dist = sqrt( (x-xp)^2 + (y-yp)^2 + (z-zp)^2 );
    %     dint = dint + h/dist;
    % end

    % dint = dint*pi*rp;
    
    % disp(['gval = ' num2str(gval)]);
    % disp(['dint = ' num2str(dint)]);
    % disp(['ratio = ' num2str(gval/dint)]);

end
