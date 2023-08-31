function [val, grad, hess] = green(k, src, targ, origin, ifdiff)
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

if nargin == 4
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

rt = repmat(targ(1,:).',1,ns);
rs = repmat(src(1,:),nt,1);
dz = repmat(src(2,:),nt,1)-repmat(targ(2,:).',1,ns);
r  = (rt + origin(1))*kabs;
dz = dz*kabs;
dr = (rs-rt)*kabs;
dout = chnk.axissymhelm2d.helm_axi(r, dr, dz, ifun, asym_tables);
int = dout.int;
intdz = dout.intdz;
intdq = dout.intdq;
intdr = dout.intdr;
intdzz = dout.intdzz;
intdrq = dout.intdrq;
intdrz = dout.intdrz;
intdqz = dout.intdqz;
for j = 1:ns
    darea = src(1,j) + origin(1);
    for i = 1:nt
        r = targ(1,i) * kabs;
        dr = (src(1,j) - targ(1,i)) * kabs;
        dz = (src(2,j) - targ(2,i)) * kabs;
        %dout = chnk.axissymhelm2d.helm_axi(r, dr, dz, ifun, asym_tables);
        if (i >1 && size(int,1)<=1)
            disp("catastrophic error")
        end
        
        vtmp(i,j) = int(i,j) * over2pi * kabs * darea;
         
        gtmp(i,j,1) = intdr(i,j) * over2pi * kabs.^2 * darea;
        gtmp(i,j,2) = intdq(i,j) * over2pi * kabs.^2 * darea;  
        gtmp(i,j,3) = -intdz(i,j) * over2pi * kabs.^2 * darea;
        
% Fix hessian scalings
% No drr, dr'r', or dr r' currently available
        %htmp(i,j,1) = 0;
        %htmp(i,j,2) = 0;
        htmp(i,j,3) = intdzz(i,j) * over2pi * kabs.^3 * darea;
        htmp(i,j,4) = intdrq(i,j) * over2pi * kabs.^3 * darea;
        htmp(i,j,5) = intdrz(i,j) * over2pi * kabs.^3 * darea;
        htmp(i,j,6) = intdqz(i,j) * over2pi * kabs.^3 * darea;
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
