function [gvals, gdzs, gdrs, gdrps] = gfuncall(zk, r, rp, dr, z, zp, dz, maxm)
%
% chnk.axissymlap2d.gfuncall evaluates a collection of axisymmetric Helmholtz
% Green's functions, defined by the expression:
%
%     gfunc(n) = \frac{rp}{4\pi} * \int_0^{2\pi} e^(i*k*|x-x'|)/|x - x'| e^(-i n t) dt 
%
% The extra factor of rp out front makes subsequent interfacing
% with RCIP slightly easier. Modes -maxm through maxm are returned, with
% gval(1) = mode -maxm, and gval(2*maxm+1) = mode maxm.
%
% The above scaling should be consistent with what is in
% chnk.axissymlap2d.gfunc, which is for merely the zero-mode
% 
    twopi = 2*pi;
    fourpi = 4*pi;
    done = 1.0;
    ima = 1i;
  
    r0 = rp;
    z0 = zp;
    rzero = sqrt(r*r + r0*r0 + dz*dz);
    alpha = 2*r*r0/rzero^2;
    x = 1/alpha;

    zkappa = zk*rzero;
    rkappa = abs(zkappa);

    athresh = 1/1.005;

  
    %!
    %! if xminus is very small, use the split method

    if (alpha > athresh)
        disp('ERROR alpha > athresh in gfuncall');
        exit
    end
    
    %!
    %! if xminus (i.e. alpha) is not small, then we can use the fft
    %! determine the number of points required
    %!
    at1 = 1/1.1;
    if (alpha >= at1)
        n = 1024;
    end
    
    if (alpha < at1)
        n = 512;
    end
    

    if (rkappa > 256) 
        nfac = log(rkappa/256)/log(2.0) + 1;
        n = 512*(2^nfac);
    end

    if (n <= maxm)
        nexp = log(maxm*done)/log(2.0);
        nexp = nexp+1;
        n = 2^nexp;
    end

    %!n = n/2

    h = 2*pi/n;


    %!
    %! construct the values
    %!
    for i = 1:n
        t = h*(i-1);
        ct = cos(t);
        rdist = rzero*sqrt(done-alpha*ct);
        zz = ima*zk*rdist;
        zzminus = zz - done;
        zzexp = exp(zz);
        denom = n* rdist^3;
        zs(i) = zzexp/(rdist)/n;
        zdr(i) = zzexp*zzminus*(r-r0*ct)/denom;
        zdr0(i) = zzexp*zzminus*(r0-r*ct)/denom;
        zdz(i) = zzexp*zzminus*(z-z0)/denom;
    end


    zs = fft(zs);
    zdr = fft(zdr);
    zdr0 = fft(zdr0);
    zdz = fft(zdz);
  
    fac = pi*rp/twopi;
    gvals = fac*zs(1:(maxm+1));
    gdrs = fac*zdr(1:(maxm+1));
    gdzs = fac*zdz(1:(maxm+1));
    gdrps = fac*zdr0(1:(maxm+1));

end
