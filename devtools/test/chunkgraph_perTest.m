%chunkgraph_perTest0();

%function chunkgraph_perTest()
%chunkgraph_perTest: 
%
% Objectives: 
% - test geometry, region detection with plot_regions,
% chunkgraph_perinregion


%housekeeping: 
clear; close all; clc; 
vrb = true;   % set false to skip figures
if vrb
    %comp domain pts:
    Nx = 150; 
    Ny = 150; 
end

%% geometry tests: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%object closed under x-tiling: 
%{
verts = [-0.5, -0.25, -0.5, 0.5, 0.25, 0.5; -4, -3.5, -3, -4, -3.5, -3]; 
edgesendverts = [6 5 1 2; 5 4 2 3]; 
merge_idx = {[1 4],[3 6]}; 
cg = chunkgraph_per(verts,edgesendverts,merge_idx); 

if vrb
    plot_geom(cg,Nx,Ny)
end

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%object closed under xy-tiling: 
%{
verts = [-0.25,-0.5,-0.5,-0.25,0.25,0.5,0.5,0.25;-5,-4.5,4.5,5,5,4.5,-4.5,-5]; 
edgesendverts = [7:-2:1;8:-2:2]; 
merge_idx = {[1 4],[2 7],[3 6],[5 8]}; 
fcurve = @(t) chnk.curves.fsine(t,0.1,pi,0); 
cg = chunkgraph_per(verts,edgesendverts,merge_idx,fcurve); 

if vrb 
    plot_geom(cg,Nx,Ny)
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%single staircase: 
%{
verts = [-0.5, -0.25, 0.25, 0.5; -0.25, 0, -0.5, -0.25];
edgesendverts = [4 3 2; 3 2 1];
merge_idx = {[1 4]};
cg = chunkgraph_per(verts,edgesendverts,merge_idx);

if vrb
    plot_geom(cg,Nx,Ny)
end

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%layered staircase
%{
verts = [-0.5, -0.25, 0.25, 0.5; -0.25, 0, -0.5, -0.25];
edgesendverts = [4 3 2; 3 2 1];
merge_idx = {[1 4]};
cg0 = chunkgraph_per(verts,edgesendverts,merge_idx);
cg  = stack_layers([cg0 + [0;-1], cg0, cg0 + [0;1]], merge_idx);

if vrb
    plot_geom(cg,Nx,Ny)
end

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
merge_idx = []; 
cg = chunkgraph_per(verts,edgesendverts,merge_idx,fchnks,cparams);

if vrb 
    plot_geom(cg,Nx,Ny); 
end
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

if vrb 
    Nx = 150; Ny = 150; 
    plot_geom(cg,Nx,Ny)
end

%}

%% chunkermat RCIP: 
%{
%QP problem on staircase: 
%
verts = [-0.5, -0.25, 0.25, 0.5; -0.25, 0, -0.5, -0.25];
edgesendverts = [4 3 2; 3 2 1];
merge_idx = {[1 4]};
cg = chunkgraph_per(verts,edgesendverts,merge_idx);

%src: 
src = []; src.r = [0;-0.5]; 

%geometry + source plot: 
figure; hold on; 
plot(cg); 
scatter(src.r(1),src.r(2),'ro','filled')

zk = 1.1; 
kern = kernel('lq','s',zk,cg.dx); 
opts = []; 

%computational domain: 
Nx = 200; Ny = 100; 
nper = 3; nshift = floor(nper/2); 
cg_comp = cg; 
dx = [cg.dx;0];
for s = -1:nshift
    if s == 0 
        continue
    end
    cg_comp = merge([cg_comp,cg + s*dx]); 
end


 
%}


%test 1: 
%{
opts = [];
opts.rcip = false;
A_no_rcip = chunkermat(cg,kern,opts);

opts.rcip = true;
A_rcip = chunkermat(cg,kern,opts);

fprintf('A_no_rcip size: %d x %d\n',size(A_no_rcip,1),size(A_no_rcip,2));
fprintf('A_rcip    size: %d x %d\n',size(A_rcip,1),size(A_rcip,2));
fprintf('relative difference = %.3e\n', ...
norm(A_rcip - A_no_rcip,'fro')/norm(A_no_rcip,'fro'));
%}

%test 2: 
%{
opts = [];
opts.rcip = true;

opts.rcip_phase = false;
[A0,~,rcip0] = chunkermat(cg,kern,opts);

opts.rcip_phase = true;
[A1,~,rcip1] = chunkermat(cg,kern,opts);

for ivert = 1:numel(rcip1)

    if isempty(rcip1{ivert})
        continue
    end

    starind = rcip1{ivert}.starind;
    pedge = rcip1{ivert}.pedge;

    ngl = cg.echnks(1).k;
    ndim = 1;  % change if vector-valued kernel

    ph = repelem(pedge(:),2*ngl*ndim);

    B_expected = (ph * (1./ph).') .* A0(starind,starind);
    B_actual   = A1(starind,starind);

    relerr = norm(B_actual - B_expected,'fro')/norm(B_expected,'fro');

    fprintf('vertex %d phase block relerr = %.3e\n',ivert,relerr);
end
%}

%test 3: 
tol = get_opt(opts,'test_tol',1e-11);

opts.rcip = true;
opts.corrections = false;
opts.nonsmoothonly = false;

fprintf('Assembling A = chunkermat(cg,kern,opts) ...\n');
[A,~,rcipsav] = chunkermat(cg,kern,opts);

rows = [];
vertex = [];
relerr = [];
abserr = [];
normref = [];
nedge_list = [];

  % ---- Local oracle comparison for each periodic representative vertex ----
    for m = 1:numel(cg.merge_idx)
        vm = cg.merge_idx{m};
        if isempty(vm), continue; end

        ivert = vm(1);  % representative/base vertex
        starind = rcipsav{ivert}.starind;
        Bgot = A(starind,starind);

        [Bref,meta] = local_periodic_rcip_reference_block(cg,kern,kappa,opts,ivert,m);

        this_abs = norm(Bgot - Bref,'fro');
        this_ref = norm(Bref,'fro');
        this_rel = this_abs/max(1,this_ref);

        fprintf('vertex %d: nedge=%d, block size=%d, relerr=%.3e, abserr=%.3e\n', ...
            ivert,meta.nedge,numel(starind),this_rel,this_abs);

        assert(this_rel < tol, ...
            'Periodic RCIP local block failed at vertex %d: relerr %.3e > %.3e.', ...
            ivert,this_rel,tol);

        rows = [rows; numel(starind)]; %#ok<AGROW>
        vertex = [vertex; ivert]; %#ok<AGROW>
        relerr = [relerr; this_rel]; %#ok<AGROW>
        abserr = [abserr; this_abs]; %#ok<AGROW>
        normref = [normref; this_ref]; %#ok<AGROW>
        nedge_list = [nedge_list; meta.nedge]; %#ok<AGROW>
    end

    results = table(vertex,nedge_list,rows,relerr,abserr,normref, ...
        'VariableNames',{'vertex','nedge','block_size','relerr','abserr','normref'});

    fprintf('All periodic RCIP local block tests passed with tol %.3e.\n',tol);

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

function plot_geom(cg,Nx,Ny)
    %basic geometry: 
    figure; hold on; 
    
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

function [kappa,kappa_p] = kappa_curve(t)
    kappa = t - 1i*sin(t); 
    kappa_p = 1 - 1i*cos(t); 
end

