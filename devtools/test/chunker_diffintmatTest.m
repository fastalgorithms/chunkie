chunker_diffintmatTest0();

function chunker_diffintmatTest0()
chnkr = chunkerfunc(@(t) [cos(t)'; 5*sin(t)']);

D = diffmat(chnkr);
C = intmat(chnkr);

dx = D * chnkr.r(1,:).';
dy = D * chnkr.r(2,:).';

% Check that (dx, dy) are unit tangent vectors
assert(vecnorm(dx.^2 + dy.^2 - 1) < 1e-10);

x = C * dx;
y = C * dy;

% Check that integral and derivative are inverse up to a scalar
assert(vecnorm(x - x(1) - chnkr.r(1,:).' + chnkr.r(1,1)) < 1e-10);
assert(vecnorm(y - y(1) - chnkr.r(2,:).' + chnkr.r(2,1)) < 1e-10);

chnkr = chunkerfunc(@(t) [cos(t)'; sin(t)']);
D = diffmat(chnkr);
test_quant = D * chnkr.r(1,:).' + chnkr.r(2,:).';
assert(vecnorm(test_quant) < 1e-10);

end