%chunkgraph_perTest0();

%function chunkgraph_perTest0()
%chunkgraph_perTest
%
% Test that chunkgraphinregion correctly labels regions for a periodic
% chunkgraph (chunkgraph_per), and that the same routine still handles an
% ordinary (non-periodic) chunkgraph.

vrb = true;   % set false to skip figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   . . . periodic chunkgraph (chunkgraph_per)

verts = [-0.5, -0.25, 0.25, 0.5; -0.25, 0, -0.5, -0.25];
edges = [4 3 2; 3 2 1];
merge_idx = {[1 4]};
cg0 = chunkgraph_per(verts,edges,merge_idx);
cg  = cg0; %stack_layers([cg0 + [0;-1], cg0], merge_idx);

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

% report (don't hard-fail on) the fraction of unclassified points; once you
% have run this and confirmed the expected value, tighten into an assert
fprintf('chunkgraph_perTest: unclassified fraction = %.4f\n', mean(isnan(ireg)));

if vrb
    figure(2); clf
    scatter(xx(:).', yy(:).', [], ireg, '.')
    axis equal; title('region ids (periodic)')

    figure(3); clf
    plot_regions(cg)
    title('plot\_regions (periodic)')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   . . . ordinary (non-periodic) chunkgraph still works through the
%         same routine (regression against the upstream code path)

verts = [[0;0],[0;1],[1;0]];
edgeendverts = [1:size(verts,2); circshift(1:size(verts,2),1)];

cgrph1 = chunkgraph(verts,edgeendverts);
cgrph2 = cgrph1 + [0;3];
cgrph  = stack_layers([cgrph1,cgrph2]);

assert(isa(cgrph,"chunkgraph") && ~isa(cgrph,"chunkgraph_per"), ...
    'stack_layers without merge_idx should return a plain chunkgraph');

x2 = linspace(-0.5,1.5,120);
y2 = linspace(-0.5,4.5,120);
[xx2,yy2] = meshgrid(x2,y2);
targs2 = []; targs2.r = [xx2(:).'; yy2(:).'];

ireg2 = chunkgraphinregion(cgrph,targs2);
lbl2  = ireg2(~isnan(ireg2));
assert(numel(ireg2) == numel(xx2), 'ireg2 should have one label per point');
assert(all(lbl2 >= 1 & lbl2 <= numel(cgrph.regions)), ...
    'region labels must lie in 1..numel(cgrph.regions)');

if vrb
    figure(4); clf
    plot(cgrph1); hold on; plot(cgrph2)
    title('two stacked (non-periodic) chunkgraphs')

    figure(5); clf
    scatter(xx2(:).', yy2(:).', [], ireg2, '.')
    axis equal; title('region ids (non-periodic)')

    figure(6); clf
    cgrph.plot_regions()
    title('plot\_regions (non-periodic)')
end

%end

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