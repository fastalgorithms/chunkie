
function gradientTest(fcn, x0, pert)

    eps = pert*randn(size(x0));

    c = 0.1;
    eps = eps/(c^2);

    [f0, g0] = fcn(x0);
    x1 = zeros(size(x0));

    fprintf("Change, result\n");
    for iter = 1:6
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
        fprintf("|eps|: %1.5e, err %1.5e\n", norm(eps,'fro'), err);
    end
end
