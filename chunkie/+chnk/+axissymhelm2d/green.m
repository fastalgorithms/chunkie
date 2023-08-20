function [val, grad, hess] = green(k, src, targ, ifdiff)
%CHNK.AXISSYMHELM2D.GREEN evaluate the helmholtz green's function
% for the given sources and targets. 
%
% Note: that the first coordinate is r, and the second z.
% The code relies on precomputed tables and hence loops are required for 
% computing various pairwise interactions.
% Finally, the code is not efficient in the sense that val, grad, hess 
% are always internally computed independent of nargout
%
% Since the kernels are not translationally invariant in r, the size
% of the gradient is 3, for d/dr, d/dr', and d/dz
%
% Similarly the hessian is of size 6 and ordered as 
%   d_{rr}, d_{r'r'}, d_{zz}, d_{rr'}, d_{rz}, d_{r'z}


zr = real(k); zi = imag(k);
if abs(zr*zi) > eps
    error('AXISSYMHELM2D.green: Axissymmetric greens function not supported for arbitrary complex', ...
           'wavenumbers, please input purely real or imaginary wave numbers\n');
end

if nargin == 3
    ifdiff = 0;
end


% Load precomputed tables
persistent asym_tables   
if isempty(asym_tables)
    asym_tables = chnk.axissymhelm2d.load_asym_tables();
end

[~, ns] = size(src);
[~, nt] = size(targ);

vtmp = zeros(nt, ns);
gtmp = zeros(nt, ns, 3);
htmp = zeros(nt, ns, 6);

% Determine whether to use 
ifun = 1;
if abs(zr) < eps
    ifun = 2;
end

if ifdiff
    if abs(zi) > eps
        error('AXISSYMHELM2D.green: Difference of greens function only supported for real wave numbers\n');
    end
    ifun = 3;
end

over2pi = 1/2/pi;
kabs = abs(k);
for j = 1:ns
    
    for i = 1:nt
        r = targ(1,i) * kabs;
        dr = (src(1,j) - targ(1,i)) * kabs;
        dz = (src(2,j) - targ(2,i)) * kabs;
        dout = chnk.axissymhelm2d.helm_axi(r, dr, dz, ifun, asym_tables);
        vtmp(i,j) = dout.int * over2pi * kabs * src(1,j);
         
        gtmp(i,j,1) = dout.intdr * over2pi * kabs.^2 * src(1,j);
        gtmp(i,j,2) = dout.intdq * over2pi * kabs.^2 * src(1,j);  
        gtmp(i,j,3) = -dout.intdz * over2pi * kabs.^2 * src(1,j);
        
% Fix hessian scalings
% No drr, dr'r', or dr r' currently available
        htmp(i,j,1) = 0;
        htmp(i,j,2) = 0;
        htmp(i,j,3) = dout.intdzz * over2pi * kabs.^3 * src(1,j);
        htmp(i,j,4) = dout.intdrq * over2pi * kabs.^3 * src(1,j);
        htmp(i,j,5) = dout.intdrz * over2pi * kabs.^3 * src(1,j);
        htmp(i,j,6) = dout.intdqz * over2pi * kabs.^3 * src(1,j);
    end
end

if nargout > 0
    val = vtmp;
end

if nargout > 1
    grad = gtmp;
end

if nargout > 2
    hess = htmp;
end
