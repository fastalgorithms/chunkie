function [gval,gdz,gdr,gdrp] = gfunc_on_axis(r,rp,dr,z,zp,dz)
    rsc = 1/(2*pi);
    n = 3;

    t = (dz.^2+dr.^2)./(2.*r.*rp);
    [q0,q1,q0d] = chnk.axissymlap2d.qleg_half(t);

    gval = rsc*(rp./r).^((n-2)/2).*q0;
    gdz  =-rsc*(rp./r).^((n-2)/2).*q0d ...
           ./(rp.*r).*dz;

    rfac = -r*(n-2)/2.*q0+(-(1+t).*r+rp).*q0d;
    gdrp  = -rsc*(rp./r).^((n-2)/2)./(rp.*r) ...
          .*rfac;
    rfac = -rp*(n-2)/2.*q0+(-(1+t).*rp+r).*q0d;
    gdr  = -rsc*(rp./r).^((n-2)/2)./(rp.*r) ...
          .*rfac;
end
