function [gvals, gdzs, gdrs, gdrps] = gfuncall(zk, r, rp, dr, z, zp, dz, maxm)
%
% chnk.axissymlap2d.gfuncall evaluates a collection of axisymmetric Laplace
% Green's functions, defined by the expression:
%
%     gfunc(n) = pi*rp * \int_0^{2\pi} e^(i*k*|x-x'|)/|x - x'| e^(-i n t) dt 
%
% The extra factor of rp (and maybe pi?) out front makes subsequent interfacing
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
  
    %r = targ(1)
    %z = targ(2)
    r0 = rp;
    z0 = zp;
    rzero = sqrt(r*r + r0*r0 + dz*dz);
    alpha = 2*r*r0/rzero^2;
    x = 1/alpha;
    % xminus = (dr*dr + dz*dz)/2/r/r0;

    zkappa = zk*rzero;
    rkappa = abs(zkappa);

    athresh = 1/1.005;

    %call prin2('x = *', x, 1)
    %call prin2('athresh = *', athresh, 1)
    %call prin2('alpha = *', alpha, 1)
    %stop
  
    %!
    %! if xminus is very small, use the split method
    %!
    %!!!!if (x .lt. 1.005d0) then
    %!!call prin2('alpha = *', alpha, 1)
    %!!call prin2('athresh = *', athresh, 1)

    if (alpha > athresh)
        %call gkmall_near(zk, src, targ, par1, maxm, vals, &
        %                 grads, grad0s)
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

    %allocate(zs(n), zdr(n), zdz(n), zdr0(n))
    %allocate(wsave(5*n))

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
        %denom = n*fourpi* rdist^3;
        denom = n* rdist^3;
        %zs(i) = zzexp/(fourpi*rdist)/n;
        zs(i) = zzexp/(rdist)/n;
        zdr(i) = zzexp*zzminus*(r-r0*ct)/denom;
        zdr0(i) = zzexp*zzminus*(r0-r*ct)/denom;
        zdz(i) = zzexp*zzminus*(z-z0)/denom;
    end

    %!
    %! take the ffts
    %!
    %call zffti(n, wsave)
    %call zfftf(n, zs, wsave)
    %call zfftf(n, zdr, wsave)
    %call zfftf(n, zdr0, wsave)
    %call zfftf(n, zdz, wsave)

    zs = fft(zs);
    zdr = fft(zdr);
    zdr0 = fft(zdr0);
    zdz = fft(zdz);
  
    %!call prin2('zs = *', zs, 2*n)
    %!call prin2('zdr = *', zdr, 2*n)
    %!call prin2('zdr0 = *', zdr0, 2*n)
    %!call prin2('zdz = *', zdz, 2*n)
    %!stop

    fac = 2*pi*rp*pi;
    gvals = fac*zs(1:(maxm+1));
    gdrs = fac*zdr(1:(maxm+1));
    gdzs = fac*zdz(1:(maxm+1));
    gdrps = fac*zdr0(1:(maxm+1));




    % check the scaling here correctly

    mode = 5;
    
    m = 1000;
    dint = 0;
    h = 2*pi/m;
    for i = 1:m
        s = (i-1)*h;
        x = r;
        y = 0;
        xp = rp*cos(s);
        yp = rp*sin(s);
        dist = sqrt( (x-xp)^2 + (y-yp)^2 + (z-zp)^2 );
        dint = dint + h*exp(ima*zk*dist)/dist*exp(-ima*mode*s);
    end
    
    dint = dint*rp*pi;
    
    disp(['gvals(mode) = ' num2str(gvals(mode+1))]);
    disp(['dint = ' num2str(dint)]);
    disp(['error = ' num2str(abs(gvals(mode+1) - dint))]);    
    disp(['ratio = ' num2str(gvals(mode+1)/dint)]);
    
end
