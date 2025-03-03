%DEMO_NEUMANN_COMBINED
% Solve the exterior Helmholtz scattering problem on a domain with corners
% and the preconditioned combined field representation
%
% Demonstrates kernel interleaving
% Demonstrates a mixture of straight and curved edges

%% Define geometry
% planewave direction
kvec = [0;8];
zk = norm(kvec);

% make vertices
narms = 5;
rots = exp(1i*2*pi*(1:narms)/narms);
verts = zeros(2,2*narms);
radi = 1;
rado = 5;
for i = 1:narms
    verts(:,2*i-1) = radi*[real(rots(i));imag(rots(i))];
    verts(:,2*i  ) = rado*[real(rots(i));imag(rots(i))];
end
nverts = size(verts,2);
edge2verts = [1:nverts;circshift(1:nverts,-1)];

% curve parameters
amp = -0.5;
frq = 5;

% define curve for each edge
fchnks = {};
for i = 1:narms
    % odd edges are straight
    fchnks{2*i-1} = [];
    % even edges are curved
    fchnks{2*i  } = @(t) sinearc(t,amp,frq);
end

cparams = [];
cparams.maxchunklen = min(4.0/zk,.5);
[cgrph] = chunkgraph(verts,edge2verts,fchnks,cparams);


% plot
figure(1);clf
plot(cgrph)
% hold on
% quiver(cgrph)
% hold off
axis equal

%% Define system

% define kernels
Sk     = kernel('helmholtz', 's', zk);
Skp    = kernel('helmholtz', 'sprime', zk);
Dk     = kernel('helmholtz', 'd', zk);

% combined field uses modified-Helmholtz kernels
Sik    = kernel('helmholtz', 's', 1i*zk);
Sikp   = kernel('helmholtz', 'sprime', 1i*zk);
Dkdiff = kernel('helmdiff', 'dprime', [zk 1i*zk]);

Z = kernel.zeros();

alpha = 1;
c1 = -1/(0.5 + 1i*alpha*0.25);
c2 = -1i*alpha/(0.5 + 1i*alpha*0.25);
c3 = -1;

% define a matrix valued kernel, which corresponds to unpacking the
% composition of operators
K = [ c1*Skp  c2*Dkdiff c2*Sikp ;
   c3*Sik  Z        Z        ;
   c3*Sikp Z        Z        ];
K = kernel(K);
Keval = c1*kernel([Sk 1i*alpha*Dk Z]);

npts = cgrph.npt;
nsys = K.opdims(1)*npts;
start = tic;
A = chunkermat(cgrph, K) + eye(nsys); 
tmat = toc(start)

%% solve system
% Set up boundary data
sources = [4;1];
strengths = 1;

rhs = zeros(nsys,1);
% get Neumann boundary data
rhs_n = -sum(kvec(:).*cgrph.n(:,:),1).*planewave(kvec(:),cgrph.r(:,:));

% only first row of K is physical
rhs(1:K.opdims(1):end) = rhs_n(:);

% solve
start = tic;
sol = gmres(A, rhs, [], 1e-13, 200);
tsolve = toc(start)

%% compute field

L = 2*rado;
x1 = linspace(-L,L,300);
[xx,yy] = meshgrid(x1,x1);
targs = [xx(:).'; yy(:).'];
ntargs = size(targs,2);

% identify points in computational domain
in = chunkerinterior(cgrph,{x1,x1});
out = ~in;

% get incoming solution
uin = nan(size(xx));
uin(out) = planewave(kvec(:),targs(:,out));

% get solution
opts = [];
opts.forcesmooth = false;
uscat = nan(size(xx));
uscat(out) = chunkerkerneval(cgrph,Keval,sol,targs(:,out),opts);

utot = uin + uscat;

%% make plots
umax = max(abs(utot(:))); 
figure(2);clf
h = pcolor(xx,yy,imag(utot)); set(h,'EdgeColor','none'); colorbar
colormap(redblue); clim([-umax,umax]);
hold on 
plot(cgrph,'k')
axis equal
title('$u^{\textrm{tot}}$','Interpreter','latex','FontSize',12)


function [r,d,d2] = sinearc(t,amp,frq)
xs = t;
ys = amp*sin(frq*t);
xp = ones(size(t));
yp = amp*frq*cos(frq*t);
xpp = zeros(size(t));
ypp = -frq*frq*amp*sin(t);

r = [(xs(:)).'; (ys(:)).'];
d = [(xp(:)).'; (yp(:)).'];
d2 = [(xpp(:)).'; (ypp(:)).'];
end
