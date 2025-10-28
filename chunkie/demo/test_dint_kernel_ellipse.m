clear


% test the volume integral of the double-layer in an ellipse
%


cparams = [];

cparams.eps = 1e-6;
cparams.nover = 0;
cparams.maxchunklen = 0.5;      



chnkr = chunkerfunc(@(t) ellipse(t,4,1), cparams);

figure(1)                                                   % plot the chunker-object 
clf
plot(chnkr, '-x')
hold on
quiver(chnkr)
axis equal


fkern = @(s,t) chnk.lap2d.kern(s,t,'d');            % double layer

opts = [];
opts.sing = 'log';

mat = chunkermat(chnkr, fkern, opts);

lhs = -(1/2).*eye(chnkr.npt) + mat   ;

nx = chnkr.n(1,:); ny = chnkr.n(2,:);

rhs = chnkr.r(1,:).^2 - chnkr.r(2,:).^2;


rhs = rhs(:);

rho = lhs\rhs;


kern_integral = @(s,t) chnk.lap2d.kern(s, t, 'dint');      % double-layer integration kernel

kern_mat = chunkermat(chnkr, kern_integral, opts);
integrand = kern_mat*rho;

numeric_result = chunkerintegral(chnkr, integrand);
analytic_result = 15*pi;
    
abserr = abs(analytic_result - numeric_result);                       % absolute error
relerr = abs((numeric_result -analytic_result)/analytic_result);      % rel error

fprintf('absolute error %5.2e\n',abserr);
fprintf('relative error %5.2e\n',relerr);