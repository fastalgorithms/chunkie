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

verts = [-0.5, -0.25, 0.25, 0.5; -0.25, 0, -0.5, -0.25]; 
edges = [4 3 2; 3 2 1]; 
merge_idx = {[1 4]}; 
cg0 = chunkgraph_per(verts,edges,merge_idx); 
cg = stack_layers([cg0 + [0;-1],cg0,cg0 + [0;1]]); 

figure
plot(cg)

%interior pts: 



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

function cgrph = stack_layers(cgrphs)
% merge a few chunkgraphs into a single chunkgraph
% assumes no vertices in common and no edges cross

nverts = 0;
verts = [];
edgesendverts = [];
fchnks = cell(0);
for i = 1:length(cgrphs)
    verts = [verts,cgrphs(i).verts];
    if class(cgrphs(i)) == "chunkgraph_per"
        edgesendverts = [edgesendverts, cgrphs(i).edgesendverts_free+nverts];
    else
        edgesendverts = [edgesendverts, cgrphs(i).edgesendverts+nverts];
    end
    nverts = nverts + size(cgrphs(i).verts,2);
    for j = 1:size(cgrphs(i).echnks,2)
        fchnks{end+1} = cgrphs(i).echnks(j);
    end
end
cgrph = chunkgraph(verts,edgesendverts,fchnks);

end