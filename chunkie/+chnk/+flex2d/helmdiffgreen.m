function [val,grad,hess,der3,der4] = helmdiffgreen(k,src,targ)
%HELMDIFFGREEN evaluate the difference of the 
% Helmholtz Green function and the Laplace Green function
% for the given sources and targets, i.e. 
%
% G(x,y) = i/4 H_0^(1)(k|x-y|) + 1/(2 pi) log(|x-y|)
%
% where H_0^(1) is the principal branch of the Hankel function
% of the first kind. This routine avoids numerical cancellation
% when |k||x-y| is small.
%
% - grad(:,:,1) has G_{x1}, grad(:,:,2) has G_{x2}
% - hess(:,:,1) has G_{x1x1}, hess(:,:,2) has G_{x1x2}, 
% hess(:,:,3) has G_{x2x2}
% - der3 has the third derivatives in the order G_{x1x1x1}, G_{x1x1x2}, 
% G_{x1x2x2}, G_{x2x2x2}
% - der4 has the fourth derivatives in the order G_{x1x1x1x1}, 
% G_{x1x1x1x2}, G_{x1x1x2x2}, G_{x1x2x2x2}, G_{x2x2x2x2}
%
% derivatives are on the *target* variables

[~,ns] = size(src);
[~,nt] = size(targ);

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);

xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);

dx = xt-xs;
dy = yt-ys;

dx2 = dx.*dx;
dx3 = dx2.*dx;
dx4 = dx3.*dx;
dy2 = dy.*dy;
dy3 = dy2.*dy;
dy4 = dy3.*dy;

r2 = dx2 + dy2;
r = sqrt(r2);
rm1 = 1./r;
rm2 = rm1.*rm1;
rm3 = rm1.*rm2;
rm4 = rm1.*rm3;
rm5 = rm1.*rm4;

% get value and r derivatives
      
[g0,g1,g2,g3,g4] = diff_h0log_and_rders(k,r);

%     evaluate potential and derivatives

if nargout > 0
    val = g0;  
end
if nargout > 1
    grad(:,:,1) = dx.*g1.*rm1;
    grad(:,:,2) = dy.*g1.*rm1;
end
if nargout > 2
    hess(:,:,1) = dx2.*g2.*rm2+g1.*(1.*rm1-dx2.*rm3);
    hess(:,:,2) = dx.*dy.*(g2.*rm2-g1.*rm3);
    hess(:,:,3) = dy2.*g2.*rm2+g1.*(1.*rm1-dy2.*rm3);
end
if nargout > 3
    der3(:,:,1) = (dx3.*g3+3*dy2.*dx.*(g2.*rm1-g1.*rm2)).*rm3;
    der3(:,:,2) = dx2.*dy.*(g3.*rm3-3*(g2.*rm4-g1.*rm5)) + ...
             dy.*(g2.*rm2-g1.*rm3);
    der3(:,:,3) = dx.*dy2.*(g3.*rm3-3*(g2.*rm4-g1.*rm5)) + ...
             dx.*(g2.*rm2-g1.*rm3);
    der3(:,:,4) = (dy3.*g3+3*dx2.*dy.*(g2.*rm1-g1.*rm2)).*rm3;
end

if nargout > 4
    der4(:,:,1) = (dx4.*(g4-6*g3.*rm1+15*(g2.*rm2-g1.*rm3))).*rm4 + ...
             (6*dx2.*(g3-3*(g2.*rm1-g1.*rm2))).*rm3 + ...
             (3*(g2-g1.*rm1)).*rm2;
    der4(:,:,2) = (dx3.*dy.*(g4-6*g3.*rm1+15*(g2.*rm2-g1.*rm3))).*rm4 + ...
             (3*dx.*dy.*(g3-3*(g2.*rm1-g1.*rm2))).*rm3;
    der4(:,:,3) = dx2.*dy2.*(g4-6*g3.*rm1+15*g2.*rm2-15*g1.*rm3).*rm4 + ...
             g3.*rm1 - 2*g2.*rm2 + 2*g1.*rm3;
    der4(:,:,4) = dx.*dy3.*(g4-6*g3.*rm1+15*(g2.*rm2-g1.*rm3)).*rm4 + ...
             3*dx.*dy.*(g3-3*(g2.*rm1-g1.*rm2)).*rm3;
    der4(:,:,5) = dy4.*(g4-6*g3.*rm1+15*(g2.*rm2-g1.*rm3)).*rm4 + ...
             6*dy2.*(g3-3*(g2.*rm1-g1.*rm2)).*rm3 + ...
             3*(g2-g1.*rm1).*rm2;
end

end

function [g0,g1,g2,g3,g4] = diff_h0log_and_rders(k,r)

g0 = zeros(size(r));
g1 = zeros(size(r));
g2 = zeros(size(r));
g3 = zeros(size(r));
g4 = zeros(size(r));

io4 = 1i*0.25;
o2p = 1/(2*pi);

isus = abs(k)*r < 1;
%isus = false(size(r));

% straightforward formulas for sufficiently large

rnot = r(~isus);

kr = k*rnot;

h0 = besselh(0,1,kr);
h1 = besselh(1,1,kr);

rm1 = 1./rnot;
rm2 = rm1.*rm1;
rm3 = rm1.*rm2;
rm4 = rm1.*rm3;

g0(~isus) = io4*h0 + o2p*log(rnot);
g1(~isus) = -k*io4*h1 + o2p*rm1;
g2(~isus) = -k*k*io4*h0 + k*io4*h1.*rm1 - o2p*rm2;
g3(~isus) = k*k*io4*h0.*rm1 + io4*k*(k*k-2*rm2).*h1 + 2*o2p*rm3;
g4(~isus) = k*io4*(3*rm2-k*k).*(2*h1.*rm1-k*h0) - 6*o2p*rm4;

% manually cancel when small

rsus = r(isus);
rm1 = 1./rsus;
rm2 = rm1.*rm1;
rm3 = rm1.*rm2;
rm4 = rm1.*rm3;

gam = 0.57721566490153286060651209;
nterms = 14;
const1 = (io4+(log(2)-gam-log(k))*o2p);

% relevant parts of hankel function represented as power series
[cf1,cf2] = chnk.flex2d.besseldiff_etc_pscoefs(nterms);
kpow = (k*k).^(1:nterms); 
cf1 = cf1(:).*kpow(:); cf2 = cf2(:).*kpow(:);

logr = log(rsus);

j0m1 = chnk.flex2d.even_pseval(cf1,rsus);
f = chnk.flex2d.even_pseval(cf2,rsus);

% differentiate power series to get derivatives
fac = 2*(1:nterms);
cf1 = cf1.*fac(:); cf2 = cf2.*fac(:);
j0m1d1 = chnk.flex2d.even_pseval(cf1,rsus).*rm1;
fd1 = chnk.flex2d.even_pseval(cf2,rsus).*rm1;
cf1 = cf1.*(fac(:)-1); cf2 = cf2.*(fac(:)-1);
j0m1d2 = chnk.flex2d.even_pseval(cf1,rsus).*rm2;

fd2 = chnk.flex2d.even_pseval(cf2,rsus).*rm2;
cf1 = cf1(:).*(fac(:)-2); cf1 = cf1(2:end);
cf2 = cf2(:).*(fac(:)-2); cf2 = cf2(2:end);
j0m1d3 = chnk.flex2d.even_pseval(cf1,rsus).*rm1;

fd3 = chnk.flex2d.even_pseval(cf2,rsus).*rm1;
fac = fac(1:end-1);
cf1 = cf1(:).*(fac(:)-1); cf2 = cf2(:).*(fac(:)-1);
j0m1d4 = chnk.flex2d.even_pseval(cf1,rsus).*rm2;
fd4 = chnk.flex2d.even_pseval(cf2,rsus).*rm2;

% combine to get derivative of i/4 H + log/(2*pi)
g0(isus) = const1*(j0m1+1) - o2p*(f  + logr.*j0m1);
g1(isus) = const1*j0m1d1 - o2p*(fd1 + logr.*j0m1d1 + j0m1.*rm1);
g2(isus) = const1*j0m1d2 - o2p*(fd2 + logr.*j0m1d2 + 2*j0m1d1.*rm1 - j0m1.*rm2);
g3(isus) = const1*j0m1d3 - o2p*(fd3 + logr.*j0m1d3 + 3*j0m1d2.*rm1 - ...
    3*j0m1d1.*rm2 + 2*j0m1.*rm3);
g4(isus) = const1*j0m1d4 - o2p*(fd4 + logr.*j0m1d4 + 4*j0m1d3.*rm1 - ...
    6*j0m1d2.*rm2 + 8*j0m1d1.*rm3 - 6*j0m1.*rm4);

end

