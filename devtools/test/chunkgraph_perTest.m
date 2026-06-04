%chunkgraph_perTest0();

%function chunkgraph_perTest0()
%chunkgraph_perTest
%
% Test that chunkgraphinregion correctly labels regions for a periodic
% chunkgraph (chunkgraph_per), and that the same routine still handles an
% ordinary (non-periodic) chunkgraph.

%housekeeping: 
clear; close all; clc; 

vrb = true;   % set false to skip figures

addpath(genpath('../../../chunkie')) % DELETE LATER

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%closed chunkgraph testing: 
%{
rng(123)
t = sort(2*pi*rand(9,1));
verts = starfish(t,5,0.3);
[~, nv] = size(verts);
edgesendverts = [1:nv; circshift(1:nv,-1)];
cpars = []; cpars.isclosed = true; 
cg = chunkgraph(verts, edgesendverts,cpars);

%pts in comp domain: 
Nx = 150; Ny = 150; 
targs = gen_comp_domain(cg,Nx,Ny); 
ireg = chunkgraphinregion(cg,targs); 

if vrb

    figure(1); hold on; 
    plot(cg); axs = axis(); hold off; 

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
    hold off; 
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%open chunkgraph testing: 
%{
verts = [-0.5, -0.25, 0.25, 0.5; -0.25, 0, -0.5, -0.25];
edgesendverts = [4 3 2; 3 2 1];
cg = chunkgraph(verts,edgesendverts);

%pts in comp domain: 
Nx = 150; Ny = 150; 
targs = gen_comp_domain(cg,Nx,Ny); 

%interior pts: 
ireg = chunkgraphinregion(cg,targs); 

if vrb
    %basic geometry: 
    figure(1); hold on; 
    plot(cg); axs = axis(); 
    title('chunkgraph\_per geometry')
    hold off; 

    %plot_regions pt id: 
    figure(2); hold on; 
    plot_regions(cg)
    axis(axs)
    title('plot\_regions() result')
    hold off; 

    %chunkgraphinregion pt id: 
    xx = targs.r(1,:); yy = targs.r(2,:); 
    nreg = size(cg.regions,2); 
    Legend = cell(1,nreg); 
    figure(3); hold on; 
    for reg = 1:nreg
        reg_idx = ireg == reg; 
        scatter(xx(reg_idx),yy(reg_idx),[],ireg(reg_idx),'.'); 
        Legend{reg} = ['region',num2str(reg)]; 
    end
    legend(Legend,'fontsize',14)
    axis(axs)
    title('chunkgraphinregion detection')
    hold off;  
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%periodic object, open in the unit cell: 
%{
verts = [-1, -0.5, -1, 1, 0.5, 1; -1, -0.5, 0, -1, -0.5, 0]; 
edgesendverts = [6 5 1 2; 5 4 2 3]; 
merge_idx = {[1 4],[3 6]}; 
cg = chunkgraph_per(verts,edgesendverts,merge_idx); 

%pts in comp domain: 
Nx = 150; Ny = 150; 
targs = gen_comp_domain(cg,Nx,Ny); 
ireg = chunkgraphinregion(cg,targs); 

if vrb

    %basic geometry: 
    figure(1); hold on; 
    plot(cg); axs = axis(); 
    title('chunkgraph\_per geometry')
    hold off; 

    %plot_regions pt id: 
    figure(2); hold on; 
    plot_regions(cg)
    axis(axs)
    title('plot\_regions() result')
    hold off; 

    %chunkgraphinregion pt id: 
    xx = targs.r(1,:); yy = targs.r(2,:); 
    nreg = size(cg.regions,2); 
    Legend = cell(1,nreg); 
    figure(3); hold on; 
    for reg = 1:nreg
        reg_idx = ireg == reg; 
        scatter(xx(reg_idx),yy(reg_idx),[],ireg(reg_idx),'.'); 
        Legend{reg} = ['region',num2str(reg)]; 
    end
    legend(Legend,'fontsize',14)
    axis(axs)
    title('chunkgraphinregion detection')
    hold off; 

end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%single staircase: 
%
verts = [-0.5, -0.25, 0.25, 0.5; -0.25, 0, -0.5, -0.25];
edgesendverts = [4 3 2; 3 2 1];
merge_idx = {[1 4]};
cg = chunkgraph_per(verts,edgesendverts,merge_idx);

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
cg  = stack_layers([cg0 + [0;-1], cg0], merge_idx);

assert(isa(cg,"chunkgraph_per"), ...
    'stack_layers with merge_idx should return a chunkgraph_per');

if vrb
    figure(1); clf
    plot(cg)
    title('chunkgraph\_per geometry')
end

% interior pts:
x1 = linspace(-0.5,0.5,150);
y1 = linspace(-1.5,1.5,150);
[xx,yy] = meshgrid(x1,y1);
targs = []; targs.r = [xx(:).'; yy(:).'];
npts = numel(xx);

ireg = chunkgraphinregion(cg,targs);
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
verts = [-0.5,-0.5,0,0,-0.5, -0.25, 0.25, 0.5; 0.5,1,1,0.5,-0.25,0,-0.5,-0.25]; 
edgesendverts = [1:4,8:-1:6;circshift(1:4,-1),7:-1:5]; 
merge_idx = {[5 8]}; 
cg = chunkgraph_per(verts,edgesendverts,merge_idx); 


%end

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

function cgrph = stack_layers(cgrphs,merge_idx)
% merge a few chunkgraphs into a single chunkgraph
% assumes no vertices in common and no edges cross

nverts = 0;
verts = [];
edgesendverts = [];
fchnks = cell(0);
verts_per = []; 
if nargin > 1
    merge_idx = repmat(merge_idx,1,numel(cgrphs)); 
end
for i = 1:length(cgrphs)
    verts = [verts,cgrphs(i).verts];
    if class(cgrphs(i)) == "chunkgraph_per"
        edgesendverts = [edgesendverts, cgrphs(i).edgesendverts_free+nverts];
        verts_per = [verts_per;cgrphs(i).vert_per]; 
        merge_idx{i} = merge_idx{i}+nverts; 
    else
        edgesendverts = [edgesendverts, cgrphs(i).edgesendverts+nverts];
    end
    nverts = nverts + size(cgrphs(i).verts,2);
    for j = 1:size(cgrphs(i).echnks,2)
        fchnks{end+1} = cgrphs(i).echnks(j);
    end
end

if nargin == 2
    cgrph = chunkgraph_per(verts,edgesendverts,merge_idx,fchnks); 

else
    cgrph = chunkgraph(verts,edgesendverts,fchnks);
end

end