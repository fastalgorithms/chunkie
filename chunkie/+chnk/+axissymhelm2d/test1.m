k = 10.1;
src = [1;2];
targ = [1.01;2];
origin = [0;0];
ifdiff = 2;
[val, grad, hess] = chnk.axissymhelm2d.green(k, src, targ, origin, ifdiff);

dx = targ-src;
dr = dx(1);
dz = dx(2);
rs = src(1);

r0   = sqrt(rs.^2+(rs+dr).^2+dz.^2);
alph = (dr.^2+dz.^2)./r0.^2;
alph = 1-alph;
frt = @(x) sqrt(1-alph*cos(x));
fex = @(x) exp(1i*k*r0*frt(x));
fdi = @(x) (fex(x)-1)./frt(x)/r0;

qint = integral(fdi,0,2*pi,'AbsTol',1E-12);
darea = rs;
qint = qint*darea/(4*pi);

dh = 1E-5;

vph = gfunc(src,targ.*[1+dh;1],k);
vmh = gfunc(src,targ.*[1-dh;1],k);
vdt1 = (vph-vmh)/(2*dh*targ(1))*darea

grad(1) - vdt1

vph = gfunc(src.*[1+dh;1],targ,k);
vmh = gfunc(src.*[1-dh;1],targ,k);
vds1 = (vph-vmh)/(2*dh*src(1))*darea

grad(2) - vds1

vph = gfunc(src,targ.*[1;1+dh],k);
vmh = gfunc(src,targ.*[1;1-dh],k);
vdt2 = (vph-vmh)/(2*dh*targ(2))*darea

grad(3) - vdt2

ds = zeros(4,4);

svec= [src;targ];
for ii=1:4
    for jj=1:4
        s = svec;
        s(ii) = s(ii)*(1+dh);
        s(jj) = s(jj)*(1+dh);
        sv = s(1:2);
        tv = s(3:4);
        vpp = gfunc(sv,tv,k);
        s = svec;
        s(ii) = s(ii)*(1-dh);
        s(jj) = s(jj)*(1+dh);
        sv = s(1:2);
        tv = s(3:4);
        vmp = gfunc(sv,tv,k); 
        s = svec;
        s(ii) = s(ii)*(1+dh);
        s(jj) = s(jj)*(1-dh);
        sv = s(1:2);
        tv = s(3:4);
        vpm = gfunc(sv,tv,k); 
        s = svec;
        s(ii) = s(ii)*(1-dh);
        s(jj) = s(jj)*(1-dh);
        sv = s(1:2);
        tv = s(3:4);
        vmm = gfunc(sv,tv,k); 
        ds(ii,jj) = (vpp-vmp+vmm-vpm)/(4*dh^2*svec(ii)*svec(jj));
    end
end


function [val] = gfunc(s_in,t_in,k)

dx = t_in-s_in;
dr = dx(1);
dz = dx(2);
rs = s_in(1);

r0   = sqrt(rs.^2+(rs+dr).^2+dz.^2);
alph = (dr.^2+dz.^2)./r0.^2;
alph = 1-alph;
frt = @(x) sqrt(1-alph*cos(x));
fex = @(x) exp(1i*k*r0*frt(x));
fdi = @(x) (fex(x)-1)./frt(x)/r0;

qint = integral(fdi,0,2*pi,'AbsTol',1E-12);
darea = 1;
val = qint*darea/(4*pi);

end