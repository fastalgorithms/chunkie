function [gval, gdz, gdr, gdrp, gdrpr, gdzz, gdrz, gdrpz] = g0funcall_vec(r, rp, dr, z, zp, dz, m)
%
% chnk.axissymlap2d.g0funcall evaluates a collection of axisymmetric Laplace
% Green's functions, defined by the expression:
%
%     gfunc(n) = pi*rp * \int_0^{2\pi} 1/|x - x'| e^(-i n t) dt
%
% Modes 0 through maxm are returned, with gval(1) = mode 0 and
% gval(maxm+1) = mode maxm.
%

    r  = reshape(r, [1,size(r,1),size(r,2)]);
    rp = reshape(rp,[1,size(rp,1),size(rp,2)]);
    dr = reshape(dr,[1,size(dr,1),size(dr,2)]);
    z  = reshape(z, [1,size(z,1),size(z,2)]);
    zp = reshape(zp,[1,size(zp,1),size(zp,2)]);
    dz = reshape(dz,[1,size(dz,1),size(dz,2)]);

    % we handle r=0 and rp=0 cases separetly 
    src_on_axis = (rp == 0);
    targ_on_axis = (r == 0) & (rp ~= 0);

    r_safe = r;
    rp_safe = rp;
    
    % dummy values whos results will be overwritten
    r_safe(r_safe == 0) = 1;
    rp_safe(rp_safe == 0) = 1;
    
    % case: rp != 0 and r != 0
    t = (dz.^2 + dr.^2)./(2.*r_safe.*rp_safe);
    chi = t + 1;

    [qm, qmd, qmdd] = chnk.axissymlap2d.qleg_half_miller_vec(t,m);

    prefac = 2*pi*sqrt(rp_safe./r_safe);

    gval = prefac.*qm;

    gdz  = prefac.*qmd ...
          ./(rp_safe.*r_safe).*dz;

    rfac = -r_safe/2.*qm + (-(1+t).*r_safe + rp_safe).*qmd;
    gdrp = prefac./(rp_safe.*r_safe).*rfac;

    rfac = -rp_safe/2.*qm + (-(1+t).*rp_safe + r_safe).*qmd;
    gdr  = prefac./(rp_safe.*r_safe).*rfac;

    rfac = 1./(rp_safe.*r_safe).*qmd ...
        + (dz./(rp_safe.*r_safe)).^2.*qmdd;
    gdzz = prefac.*rfac;

    rfac = -3./(2*r_safe.^2.*rp_safe).*qmd ...
        + (-chi./(r_safe.^2.*rp_safe) + 1./(r_safe.*rp_safe.^2)).*qmdd;
    gdrz = prefac.*dz.*rfac;

    rfac = -3./(2*rp_safe.^2.*r_safe).*qmd ...
        + (-chi./(rp_safe.^2.*r_safe) + 1./(rp_safe.*r_safe.^2)).*qmdd;
    gdrpz = prefac.*dz.*rfac;

    rfac = 1./(4*r_safe.*rp_safe).*qm ...
        + (2*chi./(rp_safe.*r_safe) ...
        - 3./(2*r_safe.^2) ...
        - 3./(2*rp_safe.^2)).*qmd ...
        + (-chi./r_safe + 1./rp_safe).*(-chi./rp_safe + 1./r_safe).*qmdd;
    gdrpr = prefac.*rfac;

    % case: rp = 0
    if any(src_on_axis(:))
        gval(:,:,src_on_axis) = 0;
        gdz(:,:,src_on_axis) = 0;
        gdr(:,:,src_on_axis) = 0;
        gdrp(:,:,src_on_axis) = 0;
        gdrpr(:,:,src_on_axis) = 0;
        gdzz(:,:,src_on_axis) = 0;
        gdrz(:,:,src_on_axis) = 0;
        gdrpz(:,:,src_on_axis) = 0;
    end
    
    % case: r = 0
    if any(targ_on_axis(:))
        gval(:,targ_on_axis) = 0;
        gdz(:,targ_on_axis) = 0;
        gdr(:,targ_on_axis) = 0;
        gdrp(:,targ_on_axis) = 0;
        gdrpr(:,targ_on_axis) = 0;
        gdzz(:,targ_on_axis) = 0;
        gdrz(:,targ_on_axis) = 0;
        gdrpz(:,targ_on_axis) = 0;

        rp0 = rp(targ_on_axis);
        dz0 = dz(targ_on_axis);

        tmp = gval(1,:,:);
        tmp(targ_on_axis) = 2*pi^2*rp0 ...
            ./ sqrt(rp0.^2 + dz0.^2);
        gval(1,:,:) = tmp;

        tmp = gdz(1,:,:);
        tmp(targ_on_axis) = -2*pi^2*rp0.*dz0 ...
            ./ sqrt(rp0.^2 + dz0.^2).^3;
        gdz(1,:,:) = tmp;

        tmp = gdrp(1,:,:);
        tmp(targ_on_axis) = -2*pi^2*rp0.^2 ...
            ./ sqrt(rp0.^2 + dz0.^2).^3;
        gdrp(1,:,:) = tmp;

        tmp = gdzz(1,:,:);
        tmp(targ_on_axis) = 2*pi^2*rp0.* ( ...
                -1 ./ sqrt(rp0.^2 + dz0.^2).^3 ...
                + 3*dz0.^2 ...
                ./ sqrt(rp0.^2 + dz0.^2).^5 ...
            );
        gdzz(1,:,:) = tmp;

        tmp = gdrpz(1,:,:);
        tmp(targ_on_axis) = 6*pi^2*rp0.^2.*dz0 ...
            ./ sqrt(rp0.^2 + dz0.^2).^5;
        gdrpz(1,:,:) = tmp;

        if m >= 1
            tmp = gdr(2,:,:);
            tmp(targ_on_axis) = pi^2*rp0.^2 ...
                ./ sqrt(rp0.^2 + dz0.^2).^3;
            gdr(2,:,:) = tmp;

            tmp = gdrz(2,:,:);
            tmp(targ_on_axis) = -3*pi^2*rp0.^2.*dz0 ...
                ./ sqrt(rp0.^2 + dz0.^2).^5;
            gdrz(2,:,:) = tmp;

            tmp = gdrpr(2,:,:);
            tmp(targ_on_axis) = pi^2*( ...
                    rp0 ./ sqrt(rp0.^2 + dz0.^2).^3 ...
                    - 3*rp0.^3 ./ sqrt(rp0.^2 + dz0.^2).^5 ...
                );
            gdrpr(2,:,:) = tmp;
        end
    end
end