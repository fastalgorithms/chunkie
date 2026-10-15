
g0funcall_vecTest0()

function g0funcall_vecTest0()
%
% verify that the modeal Green's functions for Laplace are working properly
%

h = 0.3;

r = 0.0;
rp = r + h;
dr = -h; % r - r'

z = 0.5;
zp = z-h;
dz = h; % z - z'



maxm = 10;

exact = zeros(maxm+1,1);
exact_gdz = zeros(maxm+1,1);
exact_gdr = zeros(maxm+1,1);
exact_gdrp = zeros(maxm+1,1);

exact_gdzz = zeros(maxm+1,1);
exact_gdrz = zeros(maxm+1,1);
exact_gdrpz = zeros(maxm+1,1);
exact_gdrpr = zeros(maxm+1,1);

zk = 0;
nq = 2000;
hhh = 0.0001;

% evaluate the kernel using brute force trapezoidal integration, and then
% compare with the recurrence relation code. Note: The "exact" derivatives are
% estimated using finite difference, so it is expected that the error is O(h^2)

fac = rp;

for m = 1:(maxm+1)
    dint = chnk.axissymlap2d.gkm_brute(r, rp, z, zp, zk, m-1, nq);
    exact(m) = dint*fac;

    % calculate z derivative using finite differences
    dint2 = chnk.axissymlap2d.gkm_brute(r, rp, z+hhh, zp, zk, m-1, nq);
    dint1 = chnk.axissymlap2d.gkm_brute(r, rp, z-hhh, zp, zk, m-1, nq);
    exact_gdz(m) = (dint2-dint1)/(2*hhh)*fac;
    
    % gdzz
    exact_gdzz(m) = (dint2 - 2*dint + dint1)/(hhh^2)*fac;

    % calculate r derivative using finite differences
    dint2 = chnk.axissymlap2d.gkm_brute(r+hhh, rp, z, zp, zk, m-1, nq);
    dint1 = chnk.axissymlap2d.gkm_brute(r-hhh, rp, z, zp, zk, m-1, nq);
    exact_gdr(m) = (dint2-dint1)/(2*hhh)*fac;

    % calculate rp derivative using finite differences
    dint2 = chnk.axissymlap2d.gkm_brute(r, rp+hhh, z, zp, zk, m-1, nq);
    dint1 = chnk.axissymlap2d.gkm_brute(r, rp-hhh, z, zp, zk, m-1, nq);
    exact_gdrp(m) = (dint2-dint1)/(2*hhh)*fac;
    
    % gdrz
    dpp = chnk.axissymlap2d.gkm_brute(r+hhh, rp, z+hhh, zp, zk, m-1, nq);  % +r, +z
    dpm = chnk.axissymlap2d.gkm_brute(r+hhh, rp, z-hhh, zp, zk, m-1, nq);  % +r, -z
    dmp = chnk.axissymlap2d.gkm_brute(r-hhh, rp, z+hhh, zp, zk, m-1, nq);  % -r, +z
    dmm = chnk.axissymlap2d.gkm_brute(r-hhh, rp, z-hhh, zp, zk, m-1, nq);  % -r, -z
    exact_gdrz(m) = (dpp - dpm - dmp + dmm) / (4 * hhh^2)*fac;

    % gdrpz
    dpp = chnk.axissymlap2d.gkm_brute(r, rp+hhh, z+hhh, zp, zk, m-1, nq);  % +r, +z
    dpm = chnk.axissymlap2d.gkm_brute(r, rp+hhh, z-hhh, zp, zk, m-1, nq);  % +r, -z
    dmp = chnk.axissymlap2d.gkm_brute(r, rp-hhh, z+hhh, zp, zk, m-1, nq);  % -r, +z
    dmm = chnk.axissymlap2d.gkm_brute(r, rp-hhh, z-hhh, zp, zk, m-1, nq);  % -r, -z
    exact_gdrpz(m) = (dpp - dpm - dmp + dmm) / (4 * hhh^2)*fac;

    %gdrpr
    dpp = chnk.axissymlap2d.gkm_brute(r+hhh, rp+hhh, z, zp, zk, m-1, nq);  % +r, +z
    dpm = chnk.axissymlap2d.gkm_brute(r+hhh, rp-hhh, z, zp, zk, m-1, nq);  % +r, -z
    dmp = chnk.axissymlap2d.gkm_brute(r-hhh, rp+hhh, z, zp, zk, m-1, nq);  % -r, +z
    dmm = chnk.axissymlap2d.gkm_brute(r-hhh, rp-hhh, z, zp, zk, m-1, nq);  % -r, -z
    exact_gdrpr(m) = (dpp - dpm - dmp + dmm) / (4 * hhh^2)*fac;
end


%
% compare with routine
%
[gval, gdz, gdr, gdrp, gdrpr, gdzz, gdrz, gdrpz] = chnk.axissymlap2d.g0funcall_vec(r, rp, dr, z, zp, dz, maxm);

disp(['gval error = ' num2str(norm(gval-exact))]);
disp(['gdz error = ' num2str(norm(gdz-exact_gdz))]);
disp(['gdr error = ' num2str(norm(gdr-exact_gdr))]);
disp(['gdrp error = ' num2str(norm(gdrp-exact_gdrp))]);
disp(['gdzz error = ' num2str(norm(gdzz-exact_gdzz))]);
disp(['gdrz error = ' num2str(norm(gdrz-exact_gdrz))]);
disp(['gdrpz error = ' num2str(norm(gdrpz-exact_gdrpz))]);
disp(['gdrpr error = ' num2str(norm(gdrpr-exact_gdrpr))]);

assert(norm(gval-exact) < 1e-12)
assert(max(gdz-exact_gdz) < 1e-6)
assert(max(gdr-exact_gdr) < 1e-6)
assert(max(gdrp-exact_gdrp) < 1e-6)
assert(max(gdzz-exact_gdzz) < 1e-4)
assert(max(gdrz-exact_gdrz) < 1e-4)
assert(max(gdrpz-exact_gdrpz) < 1e-4)
assert(max(gdrpr-exact_gdrpr) < 1e-4)

end