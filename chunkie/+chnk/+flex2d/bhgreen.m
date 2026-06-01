function [val,grad,hess,der3,der4,der5] = bhgreen(src,targ)
%BHGREEN evaluate the difference of the 
% Helmholtz Green function and the modified Helmholtz Green function
% for the given sources and targets, i.e. 
%
% G(x,y) = |x-y|^2 log |x-y| / 8*pi
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
%

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
      
[g0,g1,g21,g321,g4321,g54321] = r2logr_rders(r);

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

function [g0,g1,g21,g321,g4321,g54321] = r2logr_rders(r)
% g0 = g
% g1 = g'
% g21 = g'' - g'/r
% g321 = g''' - 3*g''/r + 3g'/r^2 = g''' - 3*g21/r
% g4321 = g'''' - 6*g'''/r + 15*g''/r^2 - 15*g'/r^3 
%       = g''''- 6 g321/r - 3 g21/r^2
% g54321 = g''''' - 10g''''/r + 45g'''/r^2 -105g''/r^3 + 105g'/r^4
%        = g5 - 10g4321/r - 15 g321/r^2

o8p = 1/(8*pi);

% 

r2 = r.*r;
r2d1 = 2*r;

rm1 = 1./r;
rm2 = rm1.*rm1;
rm3 = rm2.*rm1;
rm4 = rm3.*rm1;
rm5 = rm4.*rm1;

logr = log(r);

g0 = o8p*(logr.*r2);
g1 = o8p*(logr.*r2d1 + r2.*rm1);
g21 = o8p*(2*r2d1.*rm1 - 2*r2.*rm2);
g321 = o8p*(-6*r2d1.*rm2 + 8*r2.*rm3);
g4321 = o8p*(32*r2d1.*rm3 - 48*r2.*rm4);
g54321 = o8p*(- 240*r2d1.*rm4 + 384*r2.*rm5);

end

