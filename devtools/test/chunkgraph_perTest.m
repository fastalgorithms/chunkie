% chunkgraph_perTest0();
% 
% 
% function chunkgraph_perTest0()
% %chunkgraph_basicTest
% %
% % test of basic chunkgraph class constructors, methods, etc.
% %
% 



% end


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
    edgesendverts = [edgesendverts, cgrphs(i).edgesendverts+nverts];
    nverts = nverts + size(cgrphs(i).verts,2);
    for j = 1:size(cgrphs(i).echnks,2)
        fchnks{end+1} = cgrphs(i).echnks(j);
    end
end
cgrph = chunkgraph(verts,edgesendverts,fchnks);



end