%clearvars; clc;
iftorus = 0;
cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 1;
if ~iftorus % sphere
    cparams.ta = -pi/2;
    cparams.tb = pi/2;
    center = [0;2];
    cparams.ifclosed = false;
else % torus
    cparams.ta = 0;
    cparams.tb = 2*pi;
    center = [3;0];
    cparams.ifclosed = true;
end
cparams.maxchunklen = 0.5;
radius = 1;
fcurve = @(t) radius*[cos(t(:).'); sin(t(:).')];

chnkrhalf = chunkerfunc(fcurve,cparams);
chnkrhalf = chnkrhalf.move(-center); 

figure(1);clf;plot(chnkrhalf);hold on;plot(chnkrhalf,'bo');quiver(chnkrhalf,'y');

% analytic solution in n-dimension
ndim = 5;
ftrue = @(s) s.r(1,:).^2-(ndim-1)*s.r(2,:).^2; % analytical solution
fdn = @(s) -2*(ndim-1)*s.r(2,:); % (2*x1,2*x2,...,-2*(ndim-1)*z)*(0,0,1) = -2*(ndim-1)*z
fdnn = @(s) -2*(ndim-1)*ones(size(s.r(2,:)));

rhs = 2*ftrue(chnkrhalf).';
%kerns = kernel('axissymlaplace','s',ndim);
kerndp = kernel('axissymlaplace','dp',ndim);
kernd = kernel('axissymlaplace','d',ndim);
D1 = 2*kernd;
mat = chunkermat(chnkrhalf,D1);
A = mat+eye(chnkrhalf.npt);
sig = gmres(A, rhs, [], 1e-14, 200);
% place targets
x1 = linspace(center(1),center(1)+radius,200);
x2 = linspace(center(2)-radius,center(2)+radius,200);
[xx,yy] = meshgrid(x1,x2);
t = [];
t.r = [xx(:).'; yy(:).'];
in = chunkerinterior(chnkrhalf,{x1,x2});
t.r = t.r(:,in(:));
ucomp = kernd.eval(chnkrhalf,t)*(sig.*chnkrhalf.wts(:));
utrue = ftrue(t).';
err = log10(abs((ucomp-utrue)./max(abs(utrue))));

%% calculate derivative via finite difference
eps = 1e-8;
te = [];
te.r = [xx(:).'; yy(:).']+[0;1]*eps;
te.r = te.r(:,in(:));
te.n = zeros(size(t.r))+[0;1];
ucompe = kerndp.eval(chnkrhalf,te)*(sig.*chnkrhalf.wts(:));
dudn_fd = (ucompe-ucomp)./eps;
%% true derivative
dudntrue = fdn(t).';
errdudn_fd = log10(abs((dudntrue-dudn_fd)./dudntrue));

%% calculate derivative via Dp inside domain
t.n = zeros(size(t.r))+[0;1];
upcompin = kerndp.eval(chnkrhalf,t)*(sig.*chnkrhalf.wts(:));
errdudn = log10(abs((dudntrue-upcompin)./dudntrue));

%% calculate derivative via Dp on boundary
%uptrue = [2,-4]*(chnkrhalf.r(:,:).*chnkrhalf.n(:,:));
%targpt = []; targpt.r = [0;1]; targpt.n = [0;1];
%opts = []; opts.sing = 'hs';
%mat = chunkermat(chnkrhalf,Dp,opts);
%upcomp = mat*(sig.*chnkrhalf.wts(:));

%% polynomial interpolation
end1 = 1;
end2 = 1.1;
n = 40;
gridshifted = linspace(end1,end2,n);
tp = [];
tp.r = [zeros(1,n);gridshifted];
data = ftrue(tp).'; % size = (n,1)
fpoly = chebfun(data,'equi');
fprime = diff(fpoly)/(end2-end1)*2;
fprime(-1)


%% plot
tileplot = tiledlayout(1,3,'TileSpacing','compact');
plotdata1 = nan(size(xx));
plotdata1(in) = err;
plotdata2 = nan(size(xx));
plotdata2(in) = dudntrue;
plotdata3 = nan(size(xx));
plotdata3(in) = errdudn;%upcompin./dudntrue.*t.r(2,:).';


ax1 = nexttile;
plot(chnkrhalf); hold on; plot(chnkrhalf,'bo'); quiver(chnkrhalf,'r');
h = pcolor(xx,yy,reshape(plotdata1,size(xx))); set(h,'EdgeColor','none');
title('err in u'); 
clim([-10,-1]);
%clim([min(utrue),max(utrue)]);
colormap(ax1,jet(100)); colorbar; axis square;
ax2 = nexttile;
plot(chnkrhalf); hold on; plot(chnkrhalf,'bo'); quiver(chnkrhalf,'r');
h = pcolor(xx,yy,reshape(plotdata2,size(xx))); set(h,'EdgeColor','none');
title('dudntrue'); clim([min(utrue),max(utrue)]); 
clim([-18,-6]);
colormap(ax2,jet(100)); colorbar; axis square;
ax3 = nexttile;
plot(chnkrhalf); hold on; plot(chnkrhalf,'bo'); quiver(chnkrhalf,'r');
h = pcolor(xx,yy,reshape(plotdata3,size(xx))); set(h,'EdgeColor','none');
clim([-10,-1]);
%clim([1,3]);
colormap(ax3,jet(100)); colorbar; axis square;
title('upcompin');

title(tileplot,'Laplace Interior Dirichlet in n-dimension');

