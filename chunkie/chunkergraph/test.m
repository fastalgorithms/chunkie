n = 100;
[xl,wl,ul,vl]     = lege.exps(n);
[xl2,wl2,ul2,vl2] = lege.exps(n-1);

[pols,~] = lege.pols(-1,n-1);
[vfns,~] = lege.pols(xl2,n-1);

vlmod  = vfns' - ones([n-1,1])*pols';
wei_mat = 1./(1+xl2)*ones([1,n]);
vlmod(:,1) = 0;

mymat = ul2*(vlmod.*wei_mat)*ul;

x = exp(-xl);
y = exp(-(xl.^2)).*sin(xl);

xcf = mymat*x;
ycf = mymat*y;

xv = vl2*xcf;
yv = vl2*ycf;

[vfns,~] = lege.pols(xl,n-2);
xv = vfns'*xcf;
yv = vfns'*ycf;
