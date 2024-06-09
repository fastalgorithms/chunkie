%GUIDE_SIMPLEBVPS
%
% This script complements the chunkie guide section on simple boundary 
% problems.
% It shows how to discretize some common PDEs using chunkie
%

addpaths_loc();
rng(1234);

% START SPRIME
kernsp = kernel('lap','sprime');
% END SPRIME

% START LAPLACE NEUMANN
% get a chunker discretization of a starfish domain
chnkr = chunkerfunc(@(t) starfish(t));

% define a boundary condition. because of the symmetries of the 
% starfish, this function integrates to zero
pwfun = @(r) r(1,:).^2.*r(2,:);
rhs = pwfun(chnkr.r); rhs = rhs(:);

% get the kernel 
kernsp = kernel('lap','sprime');

% get a matrix discretization of the boundary integral equation 
sysmat = chunkermat(chnkr,kernsp); % just the sprime part
% add the identity term and "ones matrix"
sysmat = sysmat + 0.5*eye(chnkr.npt) + onesmat(chnkr); 

% solve the system 
sigma = gmres(sysmat,rhs);

% grid for plotting solution
x1 = linspace(-2,2,300);
[xx,yy] = meshgrid(x1,x1);
targs = [xx(:).'; yy(:).'];
in = chunkerinterior(chnkr,targs);
uu = nan(size(xx));

% need the single layer, not it's dervative, to evaluate u
kerns = kernel('lap','s');
uu(in) = chunkerkerneval(chnkr,kerns,sigma,targs(:,in));

% plot
figure(1); clf
h = pcolor(xx,yy,uu); set(h,'EdgeColor','none'); colorbar
hold on; plot(chnkr,'k'); colormap(redblue)
% END LAPLACE NEUMANN
saveas(figure(1),"guide_simplebvps_laplaceneumann.png")

%%
% START BASIC SCATTERING
% get a chunker discretization of a peanut-shaped domain
modes = [1.25,-0.25,0,0.5];
ctr = [0;0];
chnkr = chunkerfunc(@(t) chnk.curves.bymode(t,modes,ctr));

% the boundary condition is determined by an incident plane wave
kwav = 3*[-1,-5];
pwfun = @(r) exp(1i*kwav*r(:,:));
rhs = -pwfun(chnkr.r); rhs = rhs(:);

% define the CFIE kernel D - ik S
zk = norm(kwav);
coefs = [1,-1i*zk];
kerncfie = kernel('h','c',zk,coefs);

% get a matrix discretization of the boundary integral equation 
sysmat = chunkermat(chnkr,kerncfie); 
% add the identity term
sysmat = sysmat + 0.5*eye(chnkr.npt);

% solve the system 
sigma = gmres(sysmat,rhs,[],1e-10,100);

% grid for plotting solution (in exterior)
x1 = linspace(-5,5,300);
[xx,yy] = meshgrid(x1,x1);
targs = [xx(:).'; yy(:).'];
in = chunkerinterior(chnkr,{x1,x1});
uu = nan(size(xx));

% same kernel to evaluate as on boundary
uu(~in) = chunkerkerneval(chnkr,kerncfie,sigma,targs(:,~in));
uu(~in) = uu(~in) + pwfun(targs(:,~in)).';

% plot
figure(1); clf
h = pcolor(xx,yy,real(uu)); set(h,'EdgeColor','none'); colorbar
colormap(redblue); umax = max(abs(uu(:))); clim([-umax,umax]);
hold on 
plot(chnkr,'k')
% END BASIC SCATTERING
saveas(figure(1),"guide_simplebvps_basicscattering.png")

%%
% START STOKES VELOCITY
% get a chunker discretization of a peanut-shaped domain
modes = [2.5,0,0,1];
ctr = [0;0];
chnkrouter = chunkerfunc(@(t) chnk.curves.bymode(t,modes,ctr));

% inner boundaries are circles (reverse them to get orientations right)
chnkrcirc = chunkerfunc(@(t) chnk.curves.bymode(t,0.3,[0;0]));
chnkrcirc = reverse(chnkrcirc);
centers = [ [-2:2, -2:2]; [(0.7 + 0.25*(-1).^(-2:2)) , ...
    (-0.7 + 0.25*(-1).^(-2:2))]];
centers = centers + 0.1*randn(size(centers));

% make a boundary out of the outer boundary and several shifted circles
chnkrlist = [chnkrouter];
for j = 1:size(centers,2)
    chnkr1 = chnkrcirc;
    chnkr1.r(:,:) = chnkr1.r(:,:) + centers(:,j);
    chnkrlist = [chnkrlist chnkr1];
end
chnkr = merge(chnkrlist);

%%

% boundary condition specifies the velocity. two values per node
wid = 0.3; % determines width of Gaussian
f = @(r) [exp(-r(2,:).^2/(2*wid^2)); zeros(size(r(2,:)))];
rhsout = f(chnkrouter.r(:,:)); rhsout = rhsout(:);
rhs = [rhsout; zeros(2*10*chnkrcirc.npt,1)];

% define the combined layer Stokes representation kernel
mu = 1; % stokes viscosity parameter
kerndvel = kernel('stok','dvel',mu);
kernsvel = kernel('stok','svel',mu);

% get a matrix discretization of the boundary integral equation 
dmat = chunkermat(chnkr,kerndvel); 
smat = chunkermat(chnkr,kernsvel);

% add the identity term and nullspace correction
c = -1;
W = normonesmat(chnkr);
sysmat = dmat + c*smat - 0.5*eye(2*chnkr.npt) + W;

% solve the system 
sigma = gmres(sysmat,rhs,[],1e-10,100);

% grid for plotting solution (in exterior)
x1 = linspace(-3.75,3.75,300);
y1 = linspace(-2,2,300);
[xx,yy] = meshgrid(x1,y1);
targs = [xx(:).'; yy(:).'];
in = chunkerinterior(chnkr,{x1,y1});
uu = nan([2,size(xx)]);
pres = nan(size(xx));

% same kernel to evaluate as on boundary
uu(:,in) = reshape(chunkerkerneval(chnkr,kerndvel,sigma,targs(:,in)),2,nnz(in));
uu(:,in) = uu(:,in) + c*reshape(chunkerkerneval(chnkr,kernsvel,sigma,targs(:,in)),2,nnz(in));
uu(:,in) = uu(:,in) + uconst;

% can evaluate the associated pressure
kerndpres = kernel('stok','dpres',mu);
kernspres = kernel('stok','spres',mu);
opts = []; opts.eps = 1e-3;
pres(in) = chunkerkerneval(chnkr,kerndpres,sigma,targs(:,in),opts);
pres(in) = pres(in) + c*chunkerkerneval(chnkr,kernspres,sigma,targs(:,in),opts);

% plot
figure(1); clf
h = pcolor(xx,yy,pres); set(h,'EdgeColor','none'); colorbar
colormap("parula");
hold on
plot(chnkr,'k','LineWidth',2); axis equal

figure(1)
u = reshape(uu(1,:,:),size(xx)); v = reshape(uu(2,:,:),size(xx));
startt = linspace(pi-pi/12,pi+pi/12,20);
startr = 0.99*chnk.curves.bymode(startt,modes,[0;0]);
startx = startr(1,:); starty = startr(2,:);

sl = streamline(xx,yy,u,v,startx,starty);
set(sl,'LineWidth',2,'Color','w')
% END STOKES VELOCITY PROBLEM
saveas(figure(1),"guide_stokesvelocity.png")
