clearvars; close all;
addpaths_loc();

chnkr = chunkerfunc(@(t) [cos(t)'; 5*sin(t)'])

diffmat = chnkr.diffmat();
intmat = chnkr.intmat();

dx = diffmat * chnkr.r(1, 1:end)';
dy = diffmat * chnkr.r(2, 1:end)';

% Check that (dx, dy) are unit tangent vectors
assert(vecnorm(dx.^2 + dy.^2 - 1) < 1e-10);

x = intmat * dx;
y = intmat * dy;

% Check that integral and derivative are inverse up to a scalar
assert(vecnorm(x - x(1) - chnkr.r(1, 1:end)' + chnkr.r(1, 1)) < 1e-10);
assert(vecnorm(y - y(1) - chnkr.r(2, 1:end)' + chnkr.r(2, 1)) < 1e-10);

chnkr = chunkerfunc(@(t) [cos(t)'; sin(t)'])
diffmat = chnkr.diffmat();
test_quant = diffmat * chnkr.r(1, 1:end)' + chnkr.r(2, 1:end)';
assert(vecnorm(test_quant) < 1e-10);
