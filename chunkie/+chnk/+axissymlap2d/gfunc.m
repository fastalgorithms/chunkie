function [gval,gdz,gdr,gdrp] = gfunc(r,rp,dr,z,zp,dz)
    domega3 = 2*pi;
    n       = 3;
    
    t0 = (r == 0);
    s0 = (rp == 0);
    st0 = ~((~t0).*(~s0));
    dz2 = dz.^2;
    gval0 = 2*pi^2*rp./(sqrt(r.^2+rp.^2+dz2));

    gdrp0 = -2*pi*(-4*pi^2)*rp.^2./(4*pi*sqrt(rp.^2+dz.^2).^3); 
    gdz0 = -2*pi*(-4*pi^2)*rp.*dz./(4*pi*sqrt(rp.^2+dz.^2).^3); 

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

    gval(st0) = gval0(st0);
    gdr(st0) = 0;
    gdrp(st0) = gdrp0(st0);
    gdz(st0) = gdz0(st0);
end