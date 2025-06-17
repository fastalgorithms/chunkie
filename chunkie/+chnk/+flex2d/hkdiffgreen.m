function [val,grad,hess,der3,der4,der5] = hkdiffgreen(k,src,targ,ifr2logr)
%HKDIFFGREEN evaluate the difference of the 
% Helmholtz Green function and the modified Helmholtz Green function
% for the given sources and targets, i.e. 
%
% G(x,y) = i/4 H_0^(1)(k|x-y|) - 1/(2 pi) K_0(k|x-y|)
%
% or the difference of the Helmholtz and modified Helmholtz Green funcions 
% and k^2 r^2 log r/ 8 pi (a constant times the biharmonic Green function)
% i.e. 
%
% G(x,y) = i/4 H_0^(1)(k|x-y|) - 1/(2 pi) K_0(k|x-y|) + ...
%                    - k^2/(4*pi) |x-y|^2 log(|x-y|)
%
% where H_0^(1) is the principal branch of the Hankel function
% of the first kind and K_0 is the modified Bessel function of the second
% kind. This routine avoids numerical cancellation
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
%
% input:
%
% src - (2,ns) array of source locations
% targ - (2,nt) array of target locations
% k - wave number, as above
%
% optional input:
%
% ifr2logr - boolean, default: false. If true, also subtract off the 
%             k^2/(8pi) r^2 log r kernel

if nargin < 4
    ifr2logr = false;
end

r2logrfac = 1;
if ifr2logr
    r2logrfac = 0;
end

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
dx5 = dx4.*dx;

dy2 = dy.*dy;
dy3 = dy2.*dy;
dy4 = dy3.*dy;
dy5 = dy4.*dy;

r2 = dx2 + dy2;
r = sqrt(r2);
rm1 = 1./r;
rm2 = rm1.*rm1;
rm3 = rm1.*rm2;
rm4 = rm1.*rm3;
rm5 = rm1.*rm4;

% get value and r derivatives
      
[g0,g1,g21,g321,g4321,g54321] = diff_h0k0_and_rders(k,r,r2logrfac);

%     evaluate potential and derivatives

if nargout > 0
    val = g0;  
end
if nargout > 1
    grad(:,:,1) = dx.*g1.*rm1;
    grad(:,:,2) = dy.*g1.*rm1;
end
if nargout > 2
    hess(:,:,1) = dx2.*g21.*rm2+g1.*rm1;
    hess(:,:,2) = dx.*dy.*g21.*rm2;
    hess(:,:,3) = dy2.*g21.*rm2+g1.*rm1;
end
if nargout > 3
    der3(:,:,1) = dx3.*g321.*rm3+3*dx.*g21.*rm2;
    der3(:,:,2) = dx2.*dy.*g321.*rm3 + ...
             dy.*g21.*rm2;
    der3(:,:,3) = dx.*dy2.*g321.*rm3 + ...
             dx.*g21.*rm2;
    der3(:,:,4) = dy3.*g321.*rm3+3*dy.*g21.*rm2;
end

if nargout > 4
    der4(:,:,1) = dx4.*g4321.*rm4 + ...
             6*dx2.*g321.*rm3 + ...
             3*g21.*rm2;
    der4(:,:,2) = dx3.*dy.*g4321.*rm4 + ...
             3*dx.*dy.*g321.*rm3;
    der4(:,:,3) = dx2.*dy2.*g4321.*rm4 + ...
             g321.*rm1 + g21.*rm2;
    der4(:,:,4) = dx.*dy3.*g4321.*rm4 + ...
             3*dx.*dy.*g321.*rm3;
    der4(:,:,5) = dy4.*g4321.*rm4 + ...
             6*dy2.*g321.*rm3 + ...
             3*g21.*rm2;
end

if nargout > 5
    der5(:,:,1) = dx5.*g54321.*rm5 + 10*dx3.*g4321.*rm4 + 15*dx.*g321.*rm3;
    der5(:,:,2) = dy.*(dx4.*g54321.*rm5+6*dx2.*g4321.*rm4 + 3*g321.*rm3);
    der5(:,:,3) = dy2.*(dx3.*g54321.*rm5+2*dx.*g4321.*rm4) + dx.*(g4321.*rm2 + 3*g321.*rm3);
    der5(:,:,4) = dx2.*(dy3.*g54321.*rm5+2*dy.*g4321.*rm4) + dy.*(g4321.*rm2 + 3*g321.*rm3);
    der5(:,:,5) = dx.*dy4.*g54321.*rm5+6*dx.*dy2.*g4321.*rm4 + 3*dx.*g321.*rm3;
    der5(:,:,6) =  dy5.*(g54321).*rm5 + ...
        dy3.*(10*g4321).*rm4 + dy.*(15*g321).*rm3;
end

end

function [g0,g1,g21,g321,g4321,g54321] = diff_h0k0_and_rders(k,r,r2logrfac)
% g0 = g
% g1 = g'
% g21 = g'' - g'/r
% g321 = g''' - 3*g''/r + 3g'/r^2 = g''' - 3*g21/r
% g4321 = g'''' - 6*g'''/r + 15*g''/r^2 - 15*g'/r^3 
%       = g''''- 6 g321/r - 3 g21/r^2
% g54321 = g''''' - 10g''''/r + 45g'''/r^2 -105g''/r^3 + 105g'/r^4
%        = g5 - 10g4321/r - 15 g321/r^2

g0 = zeros(size(r));
g1 = zeros(size(r));
g321 = zeros(size(r));
g4321 = zeros(size(r));
g54321 = zeros(size(r));
g21 = zeros(size(r));

io4 = 1i*0.25;
o2p = 1/(2*pi);

isus = abs(k)*r < 1;
%isus = false(size(r));
%isus = true(size(r));

% straightforward formulas for sufficiently large

rnot = r(~isus);

kr = k*rnot;

h0 = besselh(0,1,kr);
h1 = besselh(1,1,kr);
h0i = besselh(0,1,1i*kr);
h1i = besselh(1,1,1i*kr);

rm1 = 1./rnot;
rm2 = rm1.*rm1;
rm3 = rm1.*rm2;
rm4 = rm1.*rm3;

r2fac = (1-r2logrfac)*k*k*0.5*o2p;
logr = log(rnot);
g0(~isus) = io4*(h0 - h0i) - r2fac*rnot.*rnot.*logr;
g1(~isus) = io4*(-k*h1 + 1i*k*h1i) - r2fac*(rnot+2*rnot.*logr);
g21(~isus) = io4*(-k*k*h0 + k*h1.*rm1 - k*k*h0i - 1i*k*h1i.*rm1) ...
    - r2fac*(3+2*logr) - g1(~isus).*rm1;
g321(~isus) = io4*(k*k*h0.*rm1 + k*(k*k-2*rm2).*h1 ...
    + k*k*h0i.*rm1 - 1i*k*(-k*k-2*rm2).*h1i) - r2fac*2*rm1 - 3*g21(~isus).*rm1;
g4321(~isus) = io4*(k*(3*rm2-k*k).*(2*h1.*rm1-k*h0) - ...
    1i*k*(3*rm2+k*k).*(2*h1i.*rm1-1i*k*h0i)) + r2fac*2*rm2 - 6*g321(~isus).*rm1 ...
    - 3*g21(~isus).*rm2;
g54321(~isus) = io4*(k*( (12*k*rm3-2*k^3*rm1).*h0  + ...
    (-24*rm4+7*k*k*rm2-k^4).*h1) - 1i*k*( (12*1i*k*rm3+1i*2*k^3*rm1).*h0i + ...
    (-24*rm4-7*k*k*rm2-k^4).*h1i)) - r2fac*4*rm3 - 10*g4321(~isus).*rm1 - ...
    15*g321(~isus).*rm2;

% manually cancel when small

rsus = r(isus);
rm1 = 1./rsus;
rm2 = rm1.*rm1;
rm3 = rm1.*rm2;
rm4 = rm1.*rm3;
rm5 = rm1.*rm4;

r2 = rsus.*rsus; 
r4 = r2.*r2;

nterms = 7;

% relevant parts of hankel function represented as power series
[cf1,cf2] = chnk.flex2d.besselkikdiff_etc_pscoefs(nterms);
cf1(1) = cf1(1)*r2logrfac;
kpow = k*k*(k.^(4*(0:(nterms-1)))); % k^2, k^6, k^10, ...
cf1 = cf1(:).*kpow(:); cf2 = cf2(:).*kpow(:);


jdiff = chnk.flex2d.pseval(cf1,r4).*r2;
f = chnk.flex2d.pseval(cf2,r4).*r2;

% differentiate power series to get derivatives
fac = 2+4*(0:(nterms-1)).';
jdiffd1 = chnk.flex2d.pseval(cf1.*fac,r4).*rsus;
fd1 = chnk.flex2d.pseval(cf2.*fac,r4).*rsus;
d = fac.*(fac-1)-fac;
jdiffd21 = chnk.flex2d.pseval(cf1.*d,r4);
fd21 = chnk.flex2d.pseval(cf2.*d,r4);
% third deriv and higher the first term is dead (this is subtle but true
% for the 21, 321, 4321, etc derivative)
cf1 = cf1(2:end); cf2 = cf2(2:end); fac = fac(2:end);
d = fac.*(fac-1).*(fac-2) - 3*(fac.*(fac-1)-fac);
jdiffd321 = chnk.flex2d.pseval(cf1.*d,r4).*rsus.*r2;
fd321 = chnk.flex2d.pseval(cf2.*d,r4).*rsus.*r2;
% g4321 = g'''' - 6*g'''/r + 15*g''/r^2 - 15*g'/r^3 
d = fac.*(fac-1).*(fac-2).*(fac-3)-6*fac.*(fac-1).*(fac-2)+15*(fac.*(fac-1)-fac);
jdiffd4321 = chnk.flex2d.pseval(cf1.*d,r4).*r2;
fd4321 = chnk.flex2d.pseval(cf2.*d,r4).*r2;
% g54321 = g''''' - 10g''''/r + 45g'''/r^2 -105g''/r^3 + 105g'/r^4
d = fac.*(fac-1).*(fac-2).*(fac-3).*(fac-4) - ...
    10*fac.*(fac-1).*(fac-2).*(fac-3) + 45*fac.*(fac-1).*(fac-2) - ...
    105*(fac.*(fac-1)-fac);
jdiffd54321 = chnk.flex2d.pseval(cf1.*d,r4).*rsus;
fd54321 = chnk.flex2d.pseval(cf2.*d,r4).*rsus;

% combine to get derivative of i/4 H - K/2pi 
gam = 0.57721566490153286060651209;
logr = log(rsus) + log(k) - log(2) + gam;
const1 = (1-r2logrfac)*o2p*(log(k)+gam-log(2))*k*k/2;

j0 = besselj(0,k*rsus);
j1 = besselj(1,k*rsus);
j2 = besselj(2,k*rsus);
j3 = besselj(3,k*rsus);
j4 = besselj(4,k*rsus);
j5 = besselj(5,k*rsus);
g0(isus) = io4*j0 - o2p*(f  + logr.*jdiff) + const1*r2;
g1(isus) = io4*(-k*j1) - o2p*(fd1 + logr.*jdiffd1 + jdiff.*rm1) + 2*const1*rsus;
g21(isus) = io4*(k*k*j2) - o2p*(fd21 + logr.*jdiffd21 + 2*jdiffd1.*rm1 - 2*jdiff.*rm2);
g321(isus) = io4*(-k^3*j3) - o2p*(fd321 + logr.*jdiffd321 + 3*jdiffd21.*rm1 - ...
    6*jdiffd1.*rm2 + 8*jdiff.*rm3);
g4321(isus) = io4*(k^4*j4)- o2p*(fd4321 + logr.*jdiffd4321 + 4*jdiffd321.*rm1 - ...
    12*jdiffd21.*rm2 + 32*jdiffd1.*rm3 - 48*jdiff.*rm4);
g54321(isus) = io4*(-k^5*j5) - o2p*(fd54321 + logr.*jdiffd54321 + ...
    5*jdiffd4321.*rm1 - 20*jdiffd321.*rm2 + 80*jdiffd21.*rm3 - 240*jdiffd1.*rm4 + 384*jdiff.*rm5);

end

