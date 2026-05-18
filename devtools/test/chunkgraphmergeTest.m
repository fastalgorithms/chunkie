chunkgraphmergeTest0();

function chunkgraphmergeTest0()
% test of chunkgraph merge, minus, and mtimes

% isosceles triangle with apex at origin, base at unit distance, half-width w
w = 0.25;
tri_verts = [0, 1, 1; 0, w, -w];
tri_edges = [1 2 3; 2 3 1];
tri = chunkgraph(tri_verts, tri_edges);

R = @(th) [cos(th) -sin(th); sin(th) cos(th)];
triA = tri;
triB = R(2*pi/3) * tri;
triC = R(4*pi/3) * tri;

vA = [1; w];
vB = R(2*pi/3)*[1; -w];
vmid = 0.5*(vA + vB);
ab = vB - vA;
perp = [-ab(2); ab(1)] / norm(ab);
vD_apex = vmid - 0.3*norm(ab)*perp;
triD = chunkgraph([vA, vB, vD_apex], [1 2 3; 2 3 1]);

cgrph = merge([triA, triB, triC, triD]);

assert(numel(cgrph.echnks) == 12);
assert(size(cgrph.verts, 2) == 8);
assert(numel(cgrph.regions) == 6);

% vector - chunkgraph uses mtimes(-1,...)
cgrph6 = [0;0] - cgrph;
assert(max(abs(cgrph6.verts(:) + cgrph.verts(:))) < 1e-13);

% plot
figure(1); clf
plot_regions(cgrph)
hold on
scatter(cgrph.verts(1,:), cgrph.verts(2,:), 60, 'r', 'filled')
axis equal
title('four merged isosceles triangles')
hold off

end
