%DEMO_MULT_CONNECT
%
% Solve the capacitor problem
% 
% Solve Laplace's equation with Dirichlet boundary condition +-1 on each
% boundary
%
% This script demonstrate an iterative solver and chunker merging


%% make geometry
cparams = [];
cparams.eps = 1.0e-6;
cparams.nover = 0;
                              
pref = []; 
pref.k = 16;
narms =5;
amp = 0.25;
scale = .3;
rad = scale*(amp+1);
ntry = 1000;

% make interior boundaries with random locations
chnkr_int = [];
L = 5;
theta = 2*pi*rand();
ctrs = L*rand()*[cos(theta);sin(theta)];
n_pts = [];
ninc = 40;
for i  = 1:ninc
    % interior boundary is rotated starfish with a random number of arms
    phi = 2*pi*rand();
    narms = randi([3,6]);
    chnkr_i = chunkerfunc(@(t) starfish(t,narms,amp,ctrs(:,i),phi,scale), ...
        cparams,pref); 
    % reverse orientation so that normals point away from computational
    % domain
    chnkr_i = reverse(chnkr_i);

    % track number of points
    n_pts(i) = chnkr_i.npt;
    chnkr_int = [chnkr_int,chnkr_i];

    % try to find location for next chunker
    for j = 1:ntry
        theta = 2*pi*rand();
        tmp = L*rand()*[cos(theta);sin(theta)];
        rmin = min(vecnorm(tmp - ctrs));
        if (rmin > rad * 3); break; end
    end
    if j == ntry; error('Could not place next boundary'); end
    ctrs = [ctrs,tmp];
end
ctrs = ctrs(:,1:end-1);
chnkr_int = merge(chnkr_int);
npt_int = chnkr_int.npt;


% make exterior boundary to be a circle containing interior boundaries
a = max(vecnorm(chnkr_int.r(:,:)))*1.1;
b = a;
chnkr_ext = chunkerfunc(@(t) ellipse(t,a,b),cparams,pref);
% refine outer curve to resolve close boundaries
chnkr_ext = refine(chnkr_ext, struct("nover",3));
n_pts(length(n_pts)+1) = chnkr_ext.npt;

chnkr = merge([chnkr_int,chnkr_ext]);

fprintf('Geometry generated\n')
figure(3); clf;
plot(chnkr)
axis equal
% hold on
% quiver(chnkr)
% hold off

%% set up system
% label indices corresponding to each boundary
idstart = [1,1+cumsum(n_pts)];

% set the Dirichlet boundary condition to be +-1 on each boundary
rhs = zeros(chnkr.npt,1);
for i = 1:length(n_pts)
    charge = (-1).^randi(2);
    rhs(idstart(i):(idstart(i+1)-1)) = charge;
end

% define representation and its gradient
% we use the combined field representation
ckern = kernel('laplace','c',[1,-1]);
ckern_grad = kernel('laplace','cgrad',[1,-1]);

tstart = tic; 
% generate quadrature corrections 
opts = [];
opts.corrections = true;
sys_cor = chunkermat(chnkr,ckern,opts);
% define system apply function
sysapply = @(sigma) -0.5*sigma + chunkermatapply(chnkr,ckern,sigma,sys_cor);

% solve system using gmres
sol = gmres(sysapply, rhs, [], 1e-10, 100); 
tsolve = toc(tstart)

%% plot electric potential

tstart = tic;
% define targets
L = max(abs(chnkr.r),[],"all");
x1 = linspace(-L,L,300);
[xx,yy] = meshgrid(x1,x1);
targs = [xx(:).'; yy(:).'];

% identify points in computational domain
in = chunkerinterior(chnkr,{x1,x1});

% evaluate potential
uu = nan(size(xx));
uu(in) = chunkerkerneval(chnkr,ckern,sol,targs(:,in));
tplot = toc(tstart)

% make plot
figure(1); clf
h = pcolor(xx,yy,uu); set(h,'EdgeColor','none'); colorbar
colormap(redblue); umax = max(abs(uu(:))); clim([-umax,umax]);
hold on 
plot(chnkr,'k')
axis equal


%% plot electric field lines

tstart = tic;
% define targets
x1 = linspace(-L,L,100);
[xx,yy] = meshgrid(x1,x1);
targs = [xx(:).'; yy(:).'];

% identify points in computational domain
in = chunkerinterior(chnkr,{x1,x1});

% evaluate electric field, the negative gradient of the potential
uu = nan(size(xx));
uugrad = nan(2,size(xx(:),1));
uugrad(:,in(:)) = -reshape(chunkerkerneval(chnkr,ckern_grad,sol,targs(:,in(:))),2,[]);
tstream = toc(tstart)

% plot fieldlines
figure(2);clf
streamslice(xx,yy,reshape(uugrad(1,:),length(x1),length(x1)),reshape(uugrad(2,:),length(x1),length(x1)),1.5)
hold on 
plot(chnkr,'k')
axis equal






