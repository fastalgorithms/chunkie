%flagrectTest
%
% tests the flagging routine based on rectangles
%

addpaths_loc();
ngrid = 100;
chnkr = chunkerfunc(@(t) starfish(t));
chnkr = refine(chnkr);
rmin = min(chnkr); rmax = max(chnkr);
dr = rmax-rmin; rmin = rmin-dr/2; rmax = rmax+dr/2;
x = linspace(rmin(1),rmax(1),ngrid);
y = linspace(rmin(2),rmax(2),ngrid);
[xx,yy] = meshgrid(x,y);

targets = [xx(:).'; yy(:).'];

sp = flagnear_rectangle(chnkr,targets);
sp2 = flagnear_rectangle_grid(chnkr,x,y);

assert(nnz(sp-sp2) == 0);

%ispec = (sp2*ones(chnkr.nch,1) > 0);

%ells = ellipses(chnkr,1.8);
%rects = bounding_rects(ells);


%clf
%plot(chnkr,'m')
%hold on
%scatter(targets(1,ispec),targets(2,ispec),'bo')
%scatter(targets(1,~ispec),targets(2,~ispec),'go')
%scatter(ells(1,:),ells(2,:),'rx')
%scatter(rects(1,:),rects(2,:),'mx')



function rects = bounding_rects(convreg)
% find minimal area bounding rectangle for each region

[dim,m,n] = size(convreg);
xreg = reshape(convreg(1,:,:),m,n);
yreg = reshape(convreg(2,:,:),m,n);

dx = diff([xreg; xreg(1,:)],1,1);
dy = diff([yreg; yreg(1,:)],1,1);
dnrm = sqrt(dx.^2+dy.^2);

dx = reshape(dx./dnrm,1,m,n);
dy = reshape(dy./dnrm,1,m,n);

d1 = [dx;dy];
d2 = [-dy;dx];

rects = zeros(2,4,n);

for i = 1:n
    d1i = d1(:,:,i);
    d2i = d2(:,:,i);
    pts = convreg(:,:,i).';
    d1c = pts*d1i;
    d2c = pts*d2i;
    d1emax = max(d1c,[],1);
    d1emin = min(d1c,[],1);
    d1eh = d1emax-d1emin;
    d2emax = max(d2c,[],1);
    d2emin = min(d2c,[],1);
    d2eh = d2emax-d2emin;
    areas = reshape(d1eh.*d2eh,m,1);
    [~,j] = min(areas);
    
    d1jmax = d1emax(j);
    d1jmin = d1emin(j);
    d2jmax = d2emax(j);
    d2jmin = d2emin(j);
    
    d1j = d1i(:,j);
    d2j = d2i(:,j);
    
    rects(:,:,i) = [d1jmax*d1j+d2jmax*d2j, d1jmin*d1j+d2jmax*d2j, ...
        d1jmin*d1j+d2jmin*d2j, d1jmax*d1j+d2jmin*d2j];
    
end
end