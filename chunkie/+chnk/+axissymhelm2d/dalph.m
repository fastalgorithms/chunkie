function [dout] = dalph(r,dr,z)

    dout = {};
    
    den = 2*r.^2+2*r.*dr+dr.^2+z.^2;
    num = 2*r.*dr+dr.^2+z.^2;
    val = num./den.^2;
    ar = -2*(r+dr).*val;
    num = -2*r.*dr-dr.^2+z.^2;
    val = num./den.^2;
    adr= -2*r.*val; 
    az = 4*(r+dr).*r.*z/den.^2;
    
    dout.ar = ar;
    dout.adr= adr;
    dout.az = az;
    
    ardr = dr.^4+4*dr.^3.*r-8*dr.*r.^3-4*r.^4-z.^4;
    ardr = 2*ardr./den.^3;

    arz = 4*(r+dr).*z.*(dr.^2+2*dr.*r-2*r.^2+z.^2)./den.^3;
    adrz= -4*r.*z.*(3*dr.^2+6*dr.*r+2*r.^2-z.^2)./den.^3;
    azz = 4*r.*(r+dr).*(dr.^2+2*dr.*r+2*r.^2-3*z.^2)./den.^3;

    dout.ardr = ardr;
    dout.arz  = arz;
    dout.adrz = adrz;
    dout.azz  = azz;
    
    r0   = sqrt(r.^2+(r+dr).^2+z.^2);
    kr   = r./r0;
    kdr  = (r+dr)./r0;
    kz   = z./r0;
    krdr = -r.*(r+dr)./r0.^3;
    krz  = -r.*z./r0.^3;
    kdrz = -(r+dr).*z./r0.^3;
    kzz  = (2*r.^2+2*r.*dr+dr.^2)./r0.^3;

    dout.kr  = kr;
    dout.kdr = kdr;
    dout.kz  = kz;
    dout.krdr= krdr;
    dout.krz = krz;
    dout.kdrz= kdrz;
    dout.kzz = kzz;
    
end

