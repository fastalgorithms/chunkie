

%
% verify that the modeal Green's functions for Laplace are working properly
%

h = 0.5;

r = 1.0;
rp = r + h;
dr = -h;

z = 0.5;
zp = z-h;
dz = h;

maxm = 30;

exact = zeros(maxm+1,1);
exact_gdz = zeros(maxm+1,1);
exact_gdr = zeros(maxm+1,1);
exact_gdrp = zeros(maxm+1,1);

zk = 0;
nq = 2000;
hhh = 0.000001

for m = 1:(maxm+1)
    dint = chnk.axissymlap2d.gkm_brute(r, rp, z, zp, zk, m-1, nq);
    exact(m) = dint*4*pi*pi*rp;

    % calculate z derivative using finite differences
    dint2 = chnk.axissymlap2d.gkm_brute(r, rp, z+hhh, zp, zk, m-1, nq);
    dint1 = chnk.axissymlap2d.gkm_brute(r, rp, z-hhh, zp, zk, m-1, nq);
    exact_gdz(m) = (dint2-dint1)/(2*hhh)*4*pi*pi*rp;

    % calculate r derivative using finite differences
    dint2 = chnk.axissymlap2d.gkm_brute(r+hhh, rp, z, zp, zk, m-1, nq);
    dint1 = chnk.axissymlap2d.gkm_brute(r-hhh, rp, z, zp, zk, m-1, nq);
    exact_gdr(m) = (dint2-dint1)/(2*hhh)*4*pi*pi*rp;

    % calculate rp derivative using finite differences
    dint2 = chnk.axissymlap2d.gkm_brute(r, rp+hhh, z, zp, zk, m-1, nq);
    dint1 = chnk.axissymlap2d.gkm_brute(r, rp-hhh, z, zp, zk, m-1, nq);
    exact_gdrp(m) = (dint2-dint1)/(2*hhh)*4*pi*pi*rp;

end

exact

exact_gdz

exact_gdr

exact_gdrp


%
% compare with routine for just zero mode
%
[gval, gdz, gdr, gdrp] = chnk.axissymlap2d.gfunc(r, rp, dr, z, zp, dz);


disp(['from 0 mode gfunc, gval = ' num2str(gval) ' error = ' num2str(abs(gval-exact(1)))]);
disp(['from 0 mode gfunc, gdz = ' num2str(gdz) ' error = ' num2str(abs(gdz-exact_gdz(1)))]);
disp(['from 0 mode gfunc, gdr = ' num2str(gdr) ' error = ' num2str(abs(gdr-exact_gdr(1)))]);
disp(['from 0 mode gfunc, gdrp = ' num2str(gdrp) ' error = ' num2str(abs(gdrp-exact_gdrp(1)))]);




% evaluate all modes now
[gvals, gdzs, gdrs, gdrps] = chnk.axissymlap2d.g0funcall(r, rp, dr, z, zp, dz, maxm);


disp('errors in all modes')

errors = abs(gvals-exact)
%errors_gdz = abs(gdzs-exact_gdz)
%errors_gdr = abs(gdrs-exact_gdr)
%errors_gdrp = abs(gdrps-exact_gdrp)
