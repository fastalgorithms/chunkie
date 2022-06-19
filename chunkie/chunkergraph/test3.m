n = 16;

[xl,wl,ul,vl] = lege.exps(n);

alpha = pi/2;
x = exp(-alpha*xl);
y = exp(-(xl.^2)).*sin(xl);

r = [x';y'];

dx = -alpha*exp(-alpha*xl);
dy = (-2*xl).*exp(-(xl.^2)).*sin(xl)+exp(-(xl.^2)).*cos(xl);

d = [dx';dy'];

ddx = alpha^2*exp(-xl);
ddy = (-2+4*xl.^2).*exp(-(xl.^2)).*sin(xl)-exp(-(xl.^2)).*sin(xl) ...
    + (-2*xl).*exp(-(xl.^2)).*cos(xl);

d2 = [ddx';ddy'];

chnkinfo.r = r;
chnkinfo.d = d;
chnkinfo.d2= d2;
[nearinfo,nearstruc] = interp2rel(chnkinfo);

nlevs = 80;
xtargs = [];
for i=1:nlevs
   xtargs = [xtargs,((xl'+1)/4+1/2)*2^(-(nlevs-i))]; 
end    

xvtrue = (exp(-alpha*((xtargs*2-1)))-exp(alpha));
xvrel  = xvtrue./xtargs;
xvexact=  -2*exp(alpha*(1-xtargs)).*sinh(alpha*xtargs);
xvrel  = xvexact./xtargs;

rnear = nearinfo.r;
xnear = rnear(1,:)';
xvv = nearstruc.eval*xnear;

xvapprx = lege.exev((xtargs*2-1),xnear).*xtargs;
xvapprxrel = xvapprx./xtargs;

figure;
plot(log(xtargs),log10(abs(xvrel-xvapprxrel)));
figure;
plot(log(xtargs),log10(abs(xvapprx-xvtrue)))


