function sing = singtype(evalfun,opdims,dim)
%CHNK.RADPP.SINGTYPE singularity type of a kernel on a curve
%
% Syntax: sing = chnk.radpp.singtype(evalfun,opdims,dim)
%
% evalfun is the kernel evaluator, evalfun(srcinfo,targinfo). The decay
% of each entry is measured on the parabola y = x^2/2 with the target at
% the origin and sources a distance rho = 1e-2 ... 1e-5 away on both
% sides. An entry that decays like rho^-q gives q = 0 for the kernels
% that are log singular or smoother, q = 1 for principal value kernels
% and q >= 2 for hypersingular kernels.
%
% see also KERNEL.RADGREEN, KERNEL.RADPP

rho = 10.^-(2:0.5:5);
nr = numel(rho);
x = [rho, -rho];
np = 2*nr;

% the points on the curve, parameterized by x
ds = sqrt(1 + x.^2);
src = struct();
src.r = zeros(dim,np);   src.r(1,:) = x;   src.r(2,:) = x.^2/2;
src.d = zeros(dim,np);   src.d(1,:) = 1;   src.d(2,:) = x;
src.d2 = zeros(dim,np);  src.d2(2,:) = 1;
src.n = zeros(dim,np);   src.n(1,:) = x./ds;   src.n(2,:) = -1./ds;
src.data = 0.7*ones(2,np);

targ = struct();
targ.r = zeros(dim,1);
targ.d = zeros(dim,1);   targ.d(1) = 1;
targ.d2 = zeros(dim,1);  targ.d2(2) = 1;
targ.n = zeros(dim,1);   targ.n(2) = -1;
targ.data = 0.7*ones(2,1);

try
    A = evalfun(src,targ);
    m = opdims(1);  n = opdims(2);
    mx = zeros(m*n,nr);
    for i = 1:m
        for j = 1:n
            v = abs(A(i,j:n:end));
            mx(i+m*(j-1),:) = max(v(1:nr),v(nr+1:end));
        end
    end
    if any(~isfinite(mx(:)))
        error('CHNK.RADPP.SINGTYPE: kernel not finite on the test curve.');
    end

    % an entry below 1e-10 of the largest at every distance is zero up to
    % rounding
    big = max(mx,[],1);
    q = 0;
    for e = 1:m*n
        if all(mx(e,:) <= 1e-10*big), continue; end
        sl = diff(log(mx(e,:) + realmin))./diff(log(rho));
        qe = max(0,ceil(-median(sl(end-2:end)) - 0.25));
        q = max(q,qe);
    end
catch ME
    warning('CHNK:radpp:singtype', ...
        'Could not determine the singularity of the kernel (%s). Using hs.', ...
        ME.message);
    q = 2;
end

if q == 0
    sing = 'log';
elseif q == 1
    sing = 'pv';
else
    sing = 'hs';
end

end
