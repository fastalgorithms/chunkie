%chunkgraph_perTest0();

function chunkgraph_perTest()
%chunkgraph_perTest
%
% Test that chunkgraphinregion correctly labels regions for a periodic
% chunkgraph (chunkgraph_per), and that the same routine still handles an
% ordinary (non-periodic) chunkgraph.

%housekeeping: 
clear; close all; clc; 
vrb = true;   % set false to skip figures

%% geometry tests: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%testing chunkgraph logic: 
%{
verts = [0,-0.5,0;0,0.5,1]; 
edgesendverts = [1:3;circshift(1:3,-1)]; 
cg = chunkgraph(verts,edgesendverts); 

figure; 
plot(cg)

figure; plot_regions(cg); 

% Generate points in the computational domain
Nx = 150; Ny = 150; 
targs = gen_comp_domain(cg, Nx, Ny); 
ireg = chunkgraphinregion(cg, targs);
xx = targs.r(1,:); yy = targs.r(2,:); 
figure; scatter(xx,yy,[],ireg,'.')
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%periodic object, open in the unit cell: 
%{
verts = [-1, -0.5, -1, 1, 0.5, 1; -1, -0.5, 0, -1, -0.5, 0]; 
edgesendverts = [6 5 1 2; 5 4 2 3]; 
merge_idx = {[1 4],[3 6]}; 
cg = chunkgraph_per(verts,edgesendverts,merge_idx); 


%basic geometry: 
figure(1); hold on; 
plot(cg); 
axs = axis(); 
title('chunkgraph\_per geometry')
hold off; 

%plot_regions verification: 
figure(3); hold on; 
plot_regions(cg)
title('plot\_regions')

%chunkgraphinregion:
Nx = 150; Ny = 150; 
targs = gen_comp_domain(cg,Nx,Ny); 
ireg = chunkgraph_perinregion(cg,targs); 
 
xx = targs.r(1,:); yy = targs.r(2,:); 
nreg = size(cg.regions,2); 
Legend = cell(1,nreg); 
figure(2); hold on; 
for reg = 1:nreg
    reg_idx = ireg == reg; 
    scatter(xx(reg_idx),yy(reg_idx),[],ireg(reg_idx),'.'); 
    Legend{reg} = ['region',num2str(reg)]; 
end
legend(Legend)
axis(axs)
title('chunkgraphinregion')
hold off; 

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%single staircase: 
%{
verts = [-0.5, -0.25, 0.25, 0.5; -0.25, 0, -0.5, -0.25];
edgesendverts = [4 3 2; 3 2 1];
merge_idx = {[1 4]};
cg = chunkgraph_per(verts,edgesendverts,merge_idx);

%basic geometry: 
figure(1); hold on; 
plot(cg); 
axs = axis(); 
title('chunkgraph\_per geometry')
hold off; 

%plot_regions: 
figure(3); hold on; 
plot_regions(cg)
title('plot\_regions')

%pts in comp domain: 
Nx = 150; Ny = 150; 
targs = gen_comp_domain(cg,Nx,Ny); 
ireg = chunkgraphinregion(cg,targs); 

%basic geometry: 
figure(1); hold on; 
plot(cg); 
axs = axis(); 
title('chunkgraph\_per geometry')
hold off; 

%chunkgraphinregion verification: 
xx = targs.r(1,:); yy = targs.r(2,:); 
nreg = size(cg.regions,2); 
Legend = cell(1,nreg); 
figure(2); hold on; 
for reg = 1:nreg
    reg_idx = ireg == reg; 
    scatter(xx(reg_idx),yy(reg_idx),[],ireg(reg_idx),'.'); 
    Legend{reg} = ['region',num2str(reg)]; 
end
legend(Legend)
axis(axs)
title('chunkgraphinregion')
hold off; 

%plot_regions verification: 
figure(3); hold on; 
plot_regions(cg)
title('plot\_regions')

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%layered staircase
%{
verts = [-0.5, -0.25, 0.25, 0.5; -0.25, 0, -0.5, -0.25];
edgesendverts = [4 3 2; 3 2 1];
merge_idx = {[1 4]};
cg0 = chunkgraph_per(verts,edgesendverts,merge_idx);
cg  = stack_layers([cg0 + [0;-1], cg0, cg0 + [0;1]], merge_idx);

%basic geometry: 
figure(1); hold on; 
plot(cg); 
axs = axis(); 
title('chunkgraph\_per geometry')
hold off; 

%plot_regions verification: 
figure(3); hold on; 
plot_regions(cg)
title('plot\_regions')

%pts in comp domain: 
Nx = 150; Ny = 150; 
targs = gen_comp_domain(cg,Nx,Ny); 
ireg = chunkgraphinregion(cg,targs); 


%chunkgraphinregion verification: 
xx = targs.r(1,:); yy = targs.r(2,:); 
nreg = size(cg.regions,2); 
Legend = cell(1,nreg); 
figure(2); hold on; 
for reg = 1:nreg
    reg_idx = ireg == reg; 
    scatter(xx(reg_idx),yy(reg_idx),[],ireg(reg_idx),'.'); 
    Legend{reg} = ['region',num2str(reg)]; 
end
legend(Legend)
axis(axs)
title('chunkgraphinregion')
hold off; 

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%closed object with user-defined periods
%{
rng(123)
t = sort(2*pi*rand(9,1));
verts = starfish(t,5,0.3);
[~, nv] = size(verts);
edgesendverts = [1:nv; circshift(1:nv,-1)];
fchnks = []; 
cparams = []; cparams.dx = 3; cparams.dy = 3; 
merge_idx = {[]}; 
cg = chunkgraph_per(verts,edgesendverts,merge_idx,fchnks,cparams);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%composite object: 
%

%closed object: 
verts = [0.4,-0.1,0.4;2,2.5,3]; 
[~, nv] = size(verts);
edgesendverts = [nv:-1:1; circshift(nv:-1:1,-1)];
fcurve = @(t) chnk.curves.fsine(t,0.1,2*pi,0); 
cpars = []; cpars.dx = []; cpars.dy = []; % dx, dy periodicities must be consistent with later geometries
merge_idx = []; 
cg1 = chunkgraph_per(verts,edgesendverts,merge_idx,fcurve,cpars); 

%object closed under x-tiling: 
verts = [-0.5, -0.25, -0.5, 0.5, 0.25, 0.5; -4, -3.5, -3, -4, -3.5, -3]; 
edgesendverts = [6 5 1 2; 5 4 2 3]; 
merge_idx = {[1 4],[3 6]}; 
cg2 = chunkgraph_per(verts,edgesendverts,merge_idx); 

%object closed under y-tiling: 
verts = [-0.2,0,0.2,-0.2,0,0.2;-5,-4.8,-5,5,4.8,5]; 
edgesendverts = [3 2 4 5;2 1 5 6]; 
merge_idx = {[1 4],[3 6]}; 
cg3 = chunkgraph_per(verts,edgesendverts,merge_idx); 

%object closed under xy-tiling: 
verts = [-0.25,-0.5,-0.5,-0.25,0.25,0.5,0.5,0.25;-5,-4.5,4.5,5,5,4.5,-4.5,-5]; 
edgesendverts = [7:-2:1;8:-2:2]; 
merge_idx = {[1 4],[2 7],[3 6],[5 8]}; 
fcurve = @(t) chnk.curves.fsine(t,0.1,pi,0); 
cg4 = chunkgraph_per(verts,edgesendverts,merge_idx,fcurve); 

%open periodic boundary: 
verts = [-0.5, -0.25, 0.25, 0.5; -0.25, 0, -0.5, -0.25];
edgesendverts = [4 3 2; 3 2 1];
merge_idx = {[1 4]};
cg5 = chunkgraph_per(verts,edgesendverts,merge_idx);
cg5 = stack_layers([cg5 + [0;-1], cg5, cg5 + [0;1], cg5 + [0;4]], merge_idx);

cg = merge([cg1,cg2,cg3,cg4,cg5]); 

%basic geometry: 
figure(1); hold on; 

subplot(1,3,1); hold on; 
plot(cg); 
axs = axis; 
title('chunkgraph\_per geometry')
hold off; 

%plot_regions: 
subplot(1,3,2); hold on; 
plot_regions(cg)
title('plot\_regions')

%chunkgraph_perinregion: 
Nx = 150; Ny = 150; 
targs = gen_comp_domain(cg,Nx,Ny); 
ireg = chunkgraph_perinregion(cg,targs); 
xx = targs.r(1,:); yy = targs.r(2,:); 
nreg = size(cg.regions,2); 
Legend = cell(1,nreg); 
subplot(1,3,3); hold on; 
for reg = 1:nreg
    reg_idx = ireg == reg; 
    scatter(xx(reg_idx),yy(reg_idx),[],ireg(reg_idx),'.'); 
    Legend{reg} = ['region',num2str(reg)]; 
end
leg = legend(Legend); leg.AutoUpdate = 'off'; 
axis(axs)
plot(cg)
title('chunkgraph\_perinregion')
hold off; 

sgtitle('region detection')
hold off; 

end

%% helpers: 
function targs = gen_comp_domain(cgrph,Nx,Ny)
    verts = cgrph.verts; 
    xmin = min(verts(1,:)); xmax = max(verts(1,:)); 
    ymin = min(verts(2,:)); ymax = max(verts(2,:)); 
    x1 = linspace(xmin-0.5,xmax+0.5,Nx);
    y1 = linspace(ymin-0.5,ymax+0.5,Ny);
    [xx,yy] = meshgrid(x1,y1);
    targs = []; targs.r = [xx(:).'; yy(:).'];
end

