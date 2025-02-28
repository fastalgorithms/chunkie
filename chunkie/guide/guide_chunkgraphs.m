%GUIDE_CHUNKGRAPHS
%
% This script complements the chunkie guide section on chunkgraph objects.
% It shows the most common methods for building and working with 
% chunkgraphs.
%

%%%%%%%%%%%%% Polygon
% START POLYGON
% chunk up random polygon

load('random_verts.mat','verts');
[~, nv] = size(verts);
edgesendverts = [1:nv; circshift(1:nv,1)];
cgrph1 = chunkgraph(verts, edgesendverts);

% plot curve, nodes, and normals

figure(1)
clf
plot(cgrph1,'b-x')
hold on
quiver(cgrph1,'r')
axis equal tight
% END POLYGON

saveas(figure(1),"guide_chunkgraphs_polygon.png");
%%
%%%%%%%%%%%%% more examples 
% START CURVED POLYGON
% Make all edges of the above chunkgraph curved


fcurve = @(t) chnk.curves.fsine(t,0.1,2*pi,0);
cgrph2 = chunkgraph(verts, edgesendverts, fcurve);


figure(2)
clf
plot(cgrph2,'b-x')
hold on
quiver(cgrph2,'r')
axis equal tight
% END CURVED POLYGON

saveas(figure(2),"guide_chunkgraphs_curved_edges.png");





%%
%%%%%%%%%%%%% more examples 
% START HAWAIIAN EARING

verts =[0;0];
funcs = cell(8,1);
funcs{1} = @(t) hawaiian(t,pi/3,0.6);
funcs{2} = @(t) hawaiian(t,pi/3,-0.6);
funcs{3} = @(t) hawaiian(t,pi/4,0.8);
funcs{4} = @(t) hawaiian(t,pi/4,-0.8);
funcs{5} = @(t) hawaiian(t,pi/6,1);
funcs{6} = @(t) hawaiian(t,pi/6,-1);

rng(1);
modes1 = 0.002*randn(11,1); modes1(1) = 0.03; 
ctr1 = [0.3;0];

modes2 = 0.1*randn(11,1); modes2(1) = 1.75; 
ctr2 = [0;0];

funcs{7} = @(t) chnk.curves.bymode(t, modes1, ctr1);
funcs{8} = @(t) chnk.curves.bymode(t, modes2, ctr2);

edgesendverts = ones(2,8);
edgesendverts(:,7:8) = nan;

cparams_closed = []; cparams_closed.ta = 0; cparams_closed.tb = 2*pi;
cparams = cell(8,1);
cparams{7} = cparams_closed;
cparams{8} = cparams_closed;
cgrph3 = chunkgraph(verts, edgesendverts, funcs, cparams);

figure(3)
clf
plot(cgrph3,'b-x')
hold on
quiver(cgrph3,'r')
axis equal tight
% END HAWAIIAN EARING

saveas(figure(3),"guide_chunkgraphs_hawaiian.png");
%%
%%%%%%%%%%%%% more examples 
% START REGION PLOT

verts = [1 0 -1 2 0 -2 0; 0 1 0 -1 2 -1 0];
edgends = [1 2 3 4 5 6 2 7; 2 3 7 5 6 4, 7, 1];
cgrph4 = chunkgraph(verts,edgends);
figure(4) 
clf
plot_regions(cgrph4);
axis equal tight

% END REGION PLOT
saveas(figure(4),"guide_chunkgraphs_regions.png");

function [r] = hawaiian(t,t0,amp)
    t = t(:).';
    rad = amp*sin(t*pi);
    ang = (t*(pi-2*t0))+t0;

    x = -rad.*sin(ang);
    y = rad.*cos(ang);

    r = [x;y];
    
end


