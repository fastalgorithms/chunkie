function gval = gkm_brute(r, rp, z, zp, zk, mode, n)

%
% This routine evaluates the modal Green's function of mode using the
% trapezoidal rule with n points:
%
%     gval = \int_0^{2\pi} e^{i k |x-x'|} / (4 pi |x-x'|) e^(-i mode t) dt 
%
% It is assumed that
%     (x,y,z) = (r,0,z) in cylindrical coordinates
%     (x',y',z') = (rp*cos(t), rp*sin(t), z') in cylindrical coordinates
%
    ima = 1j;
    
    gval = 0;
    h = 2*pi/n;
    for i = 1:n
        s = (i-1)*h;
        x = r;
        y = 0;
        xp = rp*cos(s);
        yp = rp*sin(s);
        dist = sqrt( (x-xp)^2 + (y-yp)^2 + (z-zp)^2 );
        gval = gval + exp(-ima*mode*s)*h/dist;
    end

    gval = gval/4/pi;
    

end

