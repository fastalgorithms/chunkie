function flag = flagnear_rectangle_grid(chnkr,x,y,opts)
%FLAGNEAR_RECTANGLE_GRID flag points which require special quadrature
% (or given factor of a chunklength) of any point on each chunk. on return is a sparse
% logical array. If the (i,j) entry is non zero then pts(:,i) is close
% to chunk j in chnkr.
%
% Syntax: flag = flagnear(chnkr,pts,opts)
%
% Input:
%   chnkr - chunker object describing curve
%   x - nx array of x coordinates for tensor product grid
%   y - ny array of y coordinates for tensor product grid
%
% Grid points are interpreted in the order returned by 
% MATLAB's meshgrid function [xx,yy] = meshgrid(x,y).
%
% Optional input:
%   opts - options structure
%       opts.rho = Bernstein ellipse parameter (default=1.8)
%
% Output:
%   flag - (nx*ny,chnkr.nch) sparse array. a non zero entry (i,j) means 
%       that the distance from pts(:,i) to at least one node on 
%       chunk j is less than opts.fac*length of chunk j.
%

% author: Travis Askham (askhawhat@gmail.com)

% find nearest neighbors at certain level of refinement (here chosen
% uniformly, for simplicity)

rho = 1.8;
occ = 5;

if nargin < 4
    opts = [];
end

if isfield(opts,'rho')
    rho = opts.rho;
end
if isfield(opts,'occ')
    occ = opts.occ;
end

npt = chnkr.npt;
nch = chnkr.nch;
k = chnkr.k;
dim = chnkr.dim;

% the regions needing special quadrature are images of 
% bernstein ellipses

ells = ellipses(chnkr,rho);
[rects,rectinfo] = bounding_rects(ells);

% 

[xx,yy] = meshgrid(x(:),y(:));
nx = length(x(:)); ny = length(y(:));

xmin = min(x(:)); xmax = max(x(:));
nregx = max(round(length(x)/2),1); hx = (xmax-xmin)/nregx;
xids = round((x(:)-xmin)/hx);
ymin = min(y(:)); ymax = max(y(:));
nregy = max(round(length(y)/2),1); hy = (ymax-ymin)/nregy;
yids = round((y(:)-ymin)/hy);

[xidsort,ix] = sort(xids,'ascend');
[yidsort,iy] = sort(yids,'ascend');

xifind = zeros(max(xids)+1,1);
for i = 1:length(xidsort)
    ii = xidsort(i)+1;
    xifind(ii)=xifind(ii)+1;
end
xiladr = cumsum([1;xifind(:)]);

yifind = zeros(max(yids)+1,1);
for i = 1:length(yidsort)
    ii = yidsort(i)+1;
    yifind(ii)=yifind(ii)+1;
end
yiladr = cumsum([1;yifind(:)]);

nnzero = 0;
nn = 3*nx*ny;
isp = zeros(nn,1);
jsp = zeros(nn,1);

for i = 1:nch
    rxmin = min(rects(1,:,i));
    rxmax = max(rects(1,:,i));
    rymin = min(rects(2,:,i));
    rymax = max(rects(2,:,i));
    
    ixmin = round((rxmin-xmin)/hx)+1;
    ixmax = round((rxmax-xmin)/hx)+1;
    iymin = round((rymin-ymin)/hy)+1;
    iymax = round((rymax-ymin)/hy)+1;
    
    ixrel = ix(xiladr(ixmin):(xiladr(ixmax+1)-1));
    iyrel = iy(yiladr(iymin):(yiladr(iymax+1)-1));
    
    [ixx,iyy] = meshgrid(ixrel,iyrel);
    
    xs = xx(iyrel,ixrel);
    ys = yy(iyrel,ixrel);
    
    d1 = xs*rectinfo(1,1,i)+ys*rectinfo(2,1,i);
    d2 = xs*rectinfo(1,2,i)+ys*rectinfo(2,2,i);
        
    in = and(and(and(d1 >= rectinfo(1,3,i), ...
        d1 <= rectinfo(2,3,i)), ...
        d2 >= rectinfo(1,4,i)), ...
        d2 <= rectinfo(2,4,i));
    
    it = ny*(ixx(in)-1)+iyy(in);
    nnew = nnz(in);
    js = i*ones(nnew,1);
    
    if nnew + nnzero > nn
        itemp = isp;
        jtemp = jsp;
        isp = zeros(2*nn,1);
        jsp = zeros(2*nn,1);
        isp(1:nn) = itemp;
        jsp(1:nn) = jtemp;
        nn = 2*nn;
    end
    isp(nnzero+1:nnzero+nnew) = it(:);
    jsp(nnzero+1:nnzero+nnew) = js;
    nnzero = nnew + nnzero;

end

vsp = ones(nnzero,1);
isp = isp(1:nnzero);
jsp = jsp(1:nnzero);
flag = sparse(isp,jsp,vsp,nx*ny,nch);

end

function [rects,rectinfo] = bounding_rects(convreg)
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
rectinfo = zeros(2,4,n);

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
    [~,j] = min(areas(:));
    
    d1jmax = d1emax(j);
    d1jmin = d1emin(j);
    d2jmax = d2emax(j);
    d2jmin = d2emin(j);
    
    d1j = d1i(:,j);
    d2j = d2i(:,j);
    
    rects(:,:,i) = [d1jmax*d1j+d2jmax*d2j, d1jmin*d1j+d2jmax*d2j, ...
        d1jmin*d1j+d2jmin*d2j, d1jmax*d1j+d2jmin*d2j];
    
    rectinfo(:,:,i) = [d1j,d2j,[d1jmin;d1jmax],[d2jmin;d2jmax]];
    
end

end
