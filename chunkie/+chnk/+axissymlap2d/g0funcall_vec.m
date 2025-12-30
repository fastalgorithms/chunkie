function [gval, gdz, gdr, gdrp, gdrpr, gdzz, gdrz, gdrpz] = g0funcall_vec(r, rp, dr, z, zp, dz, m)
%
% chnk.axissymlap2d.g0funcall evaluates a collection of axisymmetric Laplace
% Green's functions, defined by the expression:
%
%     gfunc(n) = pi*rp * \int_0^{2\pi} 1/|x - x'| e^(-i n t) dt
%
% it is assumed that x = (x,0,z) otherwise the integral above should pick up a
% phase factor out front of exp(i*n*phi), where phi is the azimuthal coordinate
% of x in cylindrical coordinates.
%
% The extra factor of rp (and maybe pi?) out front makes subsequent interfacing
% with RCIP slightly easier. Modes 0 through maxm are returned, with gval(1) =
% mode 0 and gval(maxm+1) = mode maxm. The function is even, so g_{-n} = g_n.
%
% The above scaling should be consistent with what is in
% chnk.axissymlap2d.gfunc, which is for merely the zero-mode
% 

    r = reshape(r,[1,size(r,1),size(r,2)]);
    rp = reshape(rp,[1,size(rp,1),size(rp,2)]);
    dr = reshape(dr,[1,size(dr,1),size(dr,2)]);
    z = reshape(z,[1,size(z,1),size(z,2)]);
    zp = reshape(zp,[1,size(zp,1),size(zp,2)]);
    dz = reshape(dz,[1,size(dz,1),size(dz,2)]);

    t = (dz.^2+dr.^2)./(2.*r.*rp);
    chi = t+1;

    [qm, qmd, qmdd] = chnk.axissymlap2d.qleg_half_miller_vec(t,m);
    
    gval = 2*pi*sqrt(rp./r).*qm;
    gdz  = 2*pi*sqrt(rp./r).*qmd ...
          ./(rp.*r).*dz;
    
    rfac = -r/2.*qm+(-(1+t).*r+rp).*qmd;
    gdrp  = 2*pi*sqrt(rp./r)./(rp.*r) ...
            .*rfac;
    
    rfac = -rp/2.*qm+(-(1+t).*rp+r).*qmd;
    gdr  = 2*pi*sqrt(rp./r)./(rp.*r) ...
           .*rfac;

    rfac = 1./(rp.*r).*qmd + (dz./(rp.*r)).^2.*qmdd;
    gdzz = 2*pi*sqrt(rp./r).*rfac;
    
    rfac = -3./(2*r.^2.*rp).*qmd + (-chi./(r.^2.*rp) + 1./(r.*rp.^2)).*qmdd;
    gdrz = 2*pi*sqrt(rp./r).*dz.*rfac;

    rfac = -3./(2*rp.^2.*r).*qmd + (-chi./(rp.^2.*r) + 1./(rp.*r.^2)).*qmdd;
    gdrpz = 2*pi*sqrt(rp./r).*dz.*rfac;

    rfac = 1./(4*r.*rp).*qm ...
        + (2*chi./(rp.*r) - 3./(2*r.^2) - 3./(2*rp.^2)).*qmd ...
        + (-chi./r+1./rp).*(-chi./rp + 1./r).*qmdd;
    gdrpr = 2*pi*sqrt(rp./r).*rfac;
end