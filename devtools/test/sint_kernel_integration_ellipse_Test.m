sint_kernel_integration_ellipse_Test0();
% test the volume integral of the single-layer in an ellipse
% 
function sint_kernel_integration_ellipse_Test0()
cparams = [];

cparams.eps = 1e-6;
cparams.nover = 0;
cparams.maxchunklen = 0.5;      



chnkr = chunkerfunc(@(t) ellipse(t,2,1), cparams);

figure(1)                                                   % plot the chunker-object 
clf
plot(chnkr, '-x')
hold on
quiver(chnkr)
axis equal


fkern = @(s,t) chnk.lap2d.kern(s,t,'sprime');

opts = [];
opts.sing = 'log';

mat = chunkermat(chnkr, fkern, opts);

lhs = (1/2).*eye(chnkr.npt) + mat + onesmat(chnkr)  ;

nx = chnkr.n(1,:); ny = chnkr.n(2,:);

rhs = 2.*chnkr.r(1,:).*nx - 2.*chnkr.r(2,:).*ny;

rhs = rhs(:);

rho = lhs\rhs;

eval_kern = @(s,t) chnk.lap2d.kern(s,t, 's');

sol = chunkerkerneval(chnkr, eval_kern, rho,[0;0]);


true_sol = 0;                    

uerr = sol - true_sol;                         % figuring out the extra constant in the interior Neumann problem



kern_integral = @(s,t) chnk.lap2d.kern(s, t, 'sint');              % single layer integration kernel

kern_mat = chunkermat(chnkr, kern_integral, opts);
integrand = kern_mat*rho;

numeric_result = chunkerintegral(chnkr, integrand);
analytic_result = 3*pi/2 + uerr*2*pi;

abserr = abs(analytic_result - numeric_result);                            % abs error
        
relerr = abs((numeric_result -analytic_result)/analytic_result);   % rel error

fprintf('absolute error %5.2e\n',abserr);
fprintf('relative error %5.2e\n',relerr);

assert(abserr<1e-9)
assert(relerr<1e-9)

end
