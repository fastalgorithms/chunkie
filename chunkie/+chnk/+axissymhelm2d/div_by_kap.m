function [dout] = div_by_kap(r,q,z,dfunc)

    r0   = sqrt(r.^2+(r+q).^2+z.^2);
    i0   =  1./r0;
    i0da =  0;
    i0dk = -1./r0.^2;
    i0daa=  0;
    i0dak=  0;
    i0dkk=  2./r0.^3;
    [dout] = chnk.axissymhelm2d.der_ak_to_grad(r,q,z,i0,i0da,i0dk,...
    i0daa,i0dak,i0dkk);

    i   = dfunc.int.*i0;
    idr = dfunc.int.*dout.intdr + dfunc.intdr.*dout.int;
    idq = dfunc.int.*dout.intdq + dfunc.intdq.*dout.int;
    idz = dfunc.int.*dout.intdz + dfunc.intdz.*dout.int;
    idqr = dfunc.intdrq.*dout.int + dfunc.int.*dout.intdrq + ...
        dfunc.intdr.*dout.intdq + dfunc.intdq.*dout.intdr;
    idrz = dfunc.intdrz.*dout.int + dfunc.int.*dout.intdrz + ...
        dfunc.intdr.*dout.intdz + dfunc.intdz.*dout.intdr;
    idqz = dfunc.intdqz.*dout.int + dfunc.int.*dout.intdqz + ...
        dfunc.intdq.*dout.intdz + dfunc.intdz.*dout.intdq;
    idzz = dfunc.intdzz.*dout.int + dfunc.int.*dout.intdzz + ...
        2*dfunc.intdz.*dout.intdz;

    dout.int   = i;
    dout.intdr = idr;
    dout.intdq = idq;
    dout.intdz = idz;
    dout.intdrq = idqr;
    dout.intdrz = idrz;
    dout.intdqz = idqz;
    dout.intdzz = idzz;
    
end

