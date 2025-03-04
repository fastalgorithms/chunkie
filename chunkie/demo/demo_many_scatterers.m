%DEMO_MANY_SCATTERERS
%
% Solve transmission scattering problems with many inclusions
%
% Demonstrate transmission problem with interleaved kernel
% Demonstrate derived quantity

% planewave direction
phi = 0;
% ambient wavenumber
zk1 = 20;
kvec = zk1*[cos(phi);sin(phi)];
% inclusion wavenumber
zk2 = 1.5*zk1;
zks = [zk1, zk2];

%% Make geometry

% geometry preferences
cparams = [];
cparams.maxchunklength = min(4.0/max(zks),0.125);
                              
pref = []; 
pref.k = 16;
narms =5;
amp = 0.25;
scale = .6;
rad = scale*(amp+1);
ntry = 1000;

% make interior boundaries with random locations
chnkr = [];
L = 3;
theta = 2*pi*rand();
ctrs = L*rand()*[cos(theta);sin(theta)];
n_pts = [];
nscat = 3;
for i  = 1:nscat
    % each boundary is rotated starfish with a random number of arms
    phi = 2*pi*rand();
    narms = randi([3,6]);
    chnkr_i = chunkerfunc(@(t) starfish(t,narms,amp,ctrs(:,i),phi,scale), ...
        cparams,pref); 
    % track number of points
    n_pts(i) = chnkr_i.npt;
    chnkr = [chnkr,chnkr_i];

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
chnkr = merge(chnkr);
npt_int = chnkr.npt;

fprintf('Geometry generated\n')
figure(3); clf;
plot(chnkr)
axis equal

%% Make system
% define system kernel

coefs = ones(2,2,2);
coefs(:,2,:) = -1;
fkern = kernel('helmdiff', 'all', zks, coefs);

% define eval kernels
fkern_eval(2,1) = kernel();
for i=1:2
    Dk = kernel('helm', 'd', zks(i));
    Sk = kernel('helm', 's', zks(i));
    fkern_eval(i) = kernel([Dk, (-1)*Sk]);
end

% define boundary data
rhs_val = -planewave(kvec,chnkr.r(:,:));
rhs_grad = 1i*sum(kvec.*chnkr.n(:,:),1).*rhs_val;

% interleave data
rhs = zeros(2*chnkr.npt,1);
rhs(1:2:end) = rhs_val;
rhs(2:2:end) = rhs_grad;

tic;
% build fast direct solver
F = chunkerflam(chnkr,fkern,1.0);
t_build_solver = toc
tic;
% solve
sol = rskelf_sv(F,rhs);
tsolve = toc


%% Compute example field

L = 1.1*max(vecnorm(chnkr.r(:,:)));
x1 = linspace(-L,L,500);
[xx,yy] = meshgrid(x1,x1);
targs = [xx(:).'; yy(:).'];
ntargs = size(targs,2);

% identify points in computational domain
in = chunkerinterior(chnkr,{x1,x1});
out = ~in;

% get incoming solution
uin = zeros(size(xx));
uin(out) = planewave(kvec(:),targs(:,out));

% get solution
tic;
uscat = zeros(size(xx));
uscat(out) = chunkerkerneval(chnkr,fkern_eval(1),sol,targs(:,out));
uscat(in) = chunkerkerneval(chnkr,fkern_eval(2),sol,targs(:,in));
tplot = toc

utot = uin + uscat;
%% make plots
umax = max(abs(utot(:))); 
figure(2);clf
h = pcolor(xx,yy,imag(utot)); set(h,'EdgeColor','none'); colorbar
colormap(redblue); clim([-umax,umax]);
hold on 
plot(chnkr,'k')
axis equal
title('$u^{\textrm{tot}}$','Interpreter','latex','FontSize',12)

%% Compute monostatic radar cross section

nphi = 500;
phis = 2*pi*(1:nphi)/nphi;

tic;
cross = zeros(nphi,1);
for i = 1:nphi
    phi = phis(i);
    kvec = zk1*[cos(phi);sin(phi)];
    
    % define boundary data
    rhs_val = -planewave(kvec,chnkr.r(:,:));
    rhs_grad = -1i*sum(kvec.*chnkr.n(:,:),1).*rhs_val;

    % interleave data
    rhs = zeros(2*chnkr.npt,1);
    rhs(1:2:end) = rhs_val;
    rhs(2:2:end) = rhs_grad;

    % solve
    sol = rskelf_sv(F,rhs);

    % extract individual densities
    ddens = sol(1:2:end);
    sdens = sol(2:2:end);
    
    % compute intermediate quantities
    d_dot_n = sum(kvec/zk1.*chnkr.n(:,:),1).';
    d_fac = exp(-1i*pi/2)*zk1;

    %compute cross section as integral against density times farfield
    %signature of kernels
    cross(i) = chnkr.wts(:).' * (sdens + ddens.*d_fac.*d_dot_n);
    cross(i) = cross(i)/4i * sqrt(zk1) * exp(-1i*pi/4) * sqrt(2/pi);
end
tcross = toc

figure(3);clf
plot(phis,abs(cross),'-','Linewidth',2)
xlabel('$\phi$','Interpreter','latex')
ylabel('abs of Monostatic cross section','Interpreter','latex')