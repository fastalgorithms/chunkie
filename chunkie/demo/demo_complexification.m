%DEMO_COMPLEXIFICATION
% Solve a Helmholtz Dirichlet scattering problem with an infinite boundary
% using the complex scaling/coordinate complexification method
%
% Demonstrates complex chunkgraphs
% Handles region detection with the complexification method

%% Problem setup
zk = 4;

% complexification parameters
eps = 1E-16;
width = 3; slope = 4;
t0 =-6; t1 = 6;
f = @(t) chnk.curves.complex_interface(t, width, slope, t0, t1);

% build geometry, two infinite rays and a triangular perturbation
xrad = -log(eps)/abs(real(zk))/slope + max([abs(t0),abs(t1)]);
endverts = real(f([-xrad;xrad]));
R = [cos(pi/4),-sin(pi/4);sin(pi/4),cos(pi/4)];
endverts(:,2) = R * endverts(:,2);

% triangle width and height
w = 1; h = -2;
verts = [endverts(:,1)-[w;0],[-w;0],[0;h],[w;0],[w;0]+endverts(:,2)];

% build solving chunkgraph
cparams = cell(1,size(verts,2)-1);
cparams{1}.ta = -xrad; cparams{1}.tb = 0;
cparams{end}.ta = 0; cparams{end}.tb = xrad;
for i = 1:length(cparams)
    cparams{i}.maxchunklen = 4/zk; 
end
cgrph = chunkgraph(verts,[1:4;2:5],{f,[],[],f},cparams); 

% build a chunkgraph with correct region detection
vertsregs = [verts, verts(:,end)-[0;10],verts(:,1)-[0;10]];
cgrphregs = chunkgraph(vertsregs,[1:7;circshift(1:7,1)]); 


%% Solve problem and evaluate fields
skern = kernel('h','s',zk);
dkern = kernel('h','d',zk);

sysmat = chunkermat(cgrph, dkern);
sysmat = sysmat - 0.5*eye(cgrph.npt);

xx = linspace(t0+6/width, t1-6/width,150);
yy = linspace(-2,4,150);
[X,Y] = meshgrid(xx,yy);
targ = []; targ.r = [X(:).';Y(:).'];

ireg = chunkgraphinregion(cgrphregs,targ);
targin = []; targin.r = targ.r(:,ireg==1);

% get incoming field
src = []; src.r = [-1;2];
rhs = -skern.eval(src,cgrph);
uin = skern.eval(src, targin);

% get scattered field, note that the fmm does not support complex points
dens = sysmat \ rhs;

opts = []; opts.accel = false;
uscat = chunkerkerneval(cgrph, dkern, dens, targin, opts);

utot = uin + uscat;

%% Plot solution
utot_fill = (NaN+1i)*zeros(size(X)); 
utot_fill(ireg==1) = utot;

figure(1);clf
h = pcolor(X,Y,imag(utot_fill)); h.EdgeColor = 'None';
hold on
plot(cgrphregs,'k-','linewidth',2)
hold off
axis equal
colorbar;
xlim([min(xx),max(xx)])
ylim([min(yy),max(yy)])
