%DEMO MIXED BOUNDARY VALUE PROBLEM
%
% This code illustrates the solution of the Laplace equation
% with mixed boundary conditions. The domain is a star-shaped
% domain where the radius has 5 Fourier modes, and the Dirichlet
% part of the boundary corresponds to y>0, and the Neumann
% part of the boundary corresponds to y<0. The boundary
% data itself is given by a Fourier series in x
% BEGIN DEMO MIXED BC
%% Define geometry
narms = 5;
amp = 0.3;

verts = [1+amp, -1+amp; 0, 0];
fchnks = @(t) starfish(t, narms, amp);
cparams = cell(2,1);
cparams{1}.ta = 0;
cparams{1}.tb = pi;
cparams{1}.maxchunklen = 0.1;
cparams{2}.ta = -pi;
cparams{2}.tb = 0;
cparams{2}.maxchunklen = 0.1;
edgesendverts = [1 2; 2 1];

cgrph = chunkgraph(verts, edgesendverts, fchnks, cparams);

%% Setup kernels
S = 2*kernel('lap', 's');
D = (-2)*kernel('lap', 'd');
Sp = 2*kernel('lap', 'sprime');
Dp = (-2)*kernel('lap', 'dprime');

K(2,2) = kernel();
K(1,1) = D;
K(1,2) = S;
K(2,1) = Dp;
K(2,2) = Sp;

Keval(1,2) = kernel();
Keval(1,1) = D;
Keval(1,2) = S;

%% Build matrix
A = chunkermat(cgrph, K);
A = A + eye(size(A,1));

%% Build test boundary data;
rhs = zeros(cgrph.npt,1);   
npt1 = cgrph.echnks(1).npt;

rng(1);
rhs(1:npt1) = real(exp(1j*40*cgrph.echnks(1).r(1,:))*exp(1j*2*pi*rand)).';
rhs(npt1+1:end) = real(exp(1j*43*cgrph.echnks(2).r(1,:))*exp(1j*2*pi*rand)).';

sig = A \ rhs;

%% Postprocess
% define targets
x1 = linspace(-1-amp,1+amp,300);
[xx,yy] = meshgrid(x1,x1);
targs = [xx(:).'; yy(:).'];

% identify points in computational domain
in = chunkerinterior(cgrph,{x1,x1});

uu = chunkerkerneval(cgrph, Keval, sig, targs);
u = nan(size(xx(:)));
u(in) = uu(in);

figure(1)
clf
plot(cgrph, 'LineWidth', 2);
hold on;
pcolor(xx, yy, reshape(u, size(xx))); shading interp
axis equal tight
% END DEMO MIXED BC
saveas(figure(1),"mixed_bc_plot.png")
