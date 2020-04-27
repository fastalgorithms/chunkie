function errs = gradient_check(fcn, x0, pert,niter,ifprint)
%GRADIENT_CHECK gradient testing function
%

    niter1 = 6;
    ifprint1 = true;
    if nargin >= 4 && ~isempty(niter)
        niter1 = niter;
    end
    if nargin >= 5 && ~isempty(ifprint)
        ifprint1 = ifprint;
    end
    
    eps = pert*randn(size(x0));

    c = 0.1;
    eps = eps/(c^2);

    [f0, g0] = fcn(x0);
    x1 = zeros(size(x0));
    
    errs = zeros(niter1,1);

    cfprintf(ifprint1,"Change, result\n");
    for iter = 1:niter1
        eps = eps*c;
        for i = 1:length(eps)
            x1(i) = x0(i) + eps(i);
        end
        [f1, g1] = fcn(x1);

        err = 0.0;
        for i = 1:length(x0)
            err = err + (g0(i) + g1(i))*eps(i);
        end
        err = 1.0 - 0.5*err/(f1 - f0);
        errs(iter) = err;
        cfprintf(ifprint1,"|eps|: %1.5e, err %1.5e\n", ...
            norm(eps,'fro'), err);
    end
end

function cfprintf(ifprint,varargin)

if (ifprint); fprintf(varargin{:}); end

end