function [valk, gradk, hessk, valik, gradik, ...
     hessik, valdiff, graddiff, hessdiff] = green_neu_all(k, src, targ, origin)
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

% Load precomputed tables
persistent asym_tables   
if isempty(asym_tables)
    asym_tables = chnk.axissymhelm2d.load_asym_tables();
end

[~, ns] = size(src);
[~, nt] = size(targ);

over2pi = 1/2/pi;
kabs = abs(k);

rt = repmat(targ(1,:).',1,ns);
rs = repmat(src(1,:),nt,1);
dz = repmat(src(2,:),nt,1)-repmat(targ(2,:).',1,ns);
r  = (rt + origin(1))*kabs;
dz = dz*kabs;
dr = (rs-rt)*kabs;
ifun = 4;
[doutk, doutik, doutdiff] = chnk.axissymhelm2d.helm_axi_all(r, dr, dz, asym_tables,ifun);

pfac = over2pi * kabs;
gfac = pfac * kabs;
hfac = gfac * kabs;
[valk, gradk, hessk] = dout_to_vgh(doutk, ns, nt, src, origin, pfac, gfac, hfac);
[valik, gradik, hessik] = dout_to_vgh(doutik, ns, nt, src, origin, pfac, gfac, hfac);
[valdiff, graddiff, hessdiff] = dout_to_vgh(doutdiff, ns, nt, src, origin, pfac, gfac, hfac);

end

function [vtmp, gtmp, htmp] = dout_to_vgh(dout, ns, nt, src, origin, pfac, gfac, hfac)
    
    vtmp = zeros(nt, ns);
    gtmp = zeros(nt, ns, 3);
    htmp = zeros(nt, ns, 6);

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
        pfacsc = pfac * darea;
        gfacsc = gfac * darea;
        hfacsc = hfac * darea;
        for i = 1:nt
            if (i >1 && size(int,1)<=1)
                disp("catastrophic error")
            end
        
            vtmp(i,j) = int(i,j) * pfacsc;
         
            gtmp(i,j,1) = intdr(i,j) * gfacsc;
            gtmp(i,j,2) = intdq(i,j) * gfacsc;  
            gtmp(i,j,3) = -intdz(i,j) * gfacsc;
        
% Fix hessian scalings
% No drr, dr'r', or dr r' currently available
        
            htmp(i,j,3) = intdzz(i,j) * hfacsc;
            htmp(i,j,4) = intdrq(i,j) * hfacsc;
            htmp(i,j,5) = intdrz(i,j) * hfacsc;
            htmp(i,j,6) = intdqz(i,j) * hfacsc;
        end
    end
end
