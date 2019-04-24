%TEST_ABSCONVGAUSS
%
% This file tests the absconvgauss routine


a = 0.75;

b = -0.75;
m = b/a;


fabs = @(x) m*abs(x);
nplot = 1000;
ntest = 10;

h = b/m/8.0;
offset = 0.0;

xplot = linspace(-a,a,nplot);

figure(1)
clf
plot(xplot,fabs(xplot),'r')
hold on
plot(xplot,chnk.spcl.absconvgauss(xplot,m,offset,h));

figure(2)
clf
ff = chnk.spcl.absconvgauss(xplot,m,offset,h);
df = abs(ff-fabs(xplot));
semilogy(df);

xtest = linspace(-a,a,ntest);
pert = 0.1*b/m;

for i = 1:ntest
    x0 = xtest(i);
    gradientTest(@(x) chnk.spcl.absconvgauss(x,m,offset,h),x0,pert);
    gradientTest(@(x) acgder(x,m,offset,h),x0,pert);
end

function [d,d2] = acgder(x,a,b,h)
[~,d,d2] = chnk.spcl.absconvgauss(x,a,b,h);
end
