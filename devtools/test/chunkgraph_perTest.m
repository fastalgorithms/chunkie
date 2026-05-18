% chunkgraph_perTest0();
% 
% 
% function chunkgraph_perTest0()
% %chunkgraph_basicTest
% %
% % test of basic chunkgraph class constructors, methods, etc.
% %
% 

%housekeeping: 
clear all; close all; clc; 

%filepath: 
addpath(genpath('../../../../chunkie'))

verts = [-0.5, -0.25, 0.25, 0.5; -0.25, 0, -0.5, -0.25]; 
edges = [4 3 2; 3 2 1]; 
merge_idx = {[1 4]}; 
cg0 = chunkgraph_per(verts,edges,merge_idx); 
cg = cg0; 
cg = stack_layers([cg0 + [0;-1],cg0],merge_idx); 

figure
plot(cg)

%interior pts: 
x1 = linspace(-0.5,0.5,150); 
y1 = linspace(-1.5,1.5,150); 
[xx,yy] = meshgrid(x1,y1); 
targs = []; targs.r = [xx(:).';yy(:).']; 

ireg = chunkgraphinregion(cg,targs); 
figure
scatter(xx(:).',yy(:).',[],ireg,'.')

figure
plot_regions(cg)

verts = [[0;0],[0;1],[1;0]];

edgeendverts = [1:size(verts,2);circshift(1:size(verts,2),1)];

cgrph1 = chunkgraph(verts,edgeendverts);
cgrph2 = cgrph1 + [0;3];

figure(1);clf
plot(cgrph1)
hold on
plot(cgrph2)


cgrph = stack_layers([cgrph1,cgrph2]);

figure(2)
cgrph.plot_regions()

function cgrph = stack_layers(cgrphs,merge_idx)
% merge a few chunkgraphs into a single chunkgraph
% assumes no vertices in common and no edges cross

nverts = 0;
verts = [];
edgesendverts = [];
fchnks = cell(0);
verts_per = []; 
merge_idx = repmat(merge_idx,1,numel(cgrphs)); 
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