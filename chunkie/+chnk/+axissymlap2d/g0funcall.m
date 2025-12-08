function [gvals, gdzs, gdrs, gdrps] = g0funcall(r, rp, dr, z, zp, dz, maxm)

%
% chnk.axissymlap2d.g0funcall evaluates a collection of axisymmetric Laplace
% Green's functions, defined by the expression:
%
%     gfunc(n) = pi*rp * \int_0^{2\pi} 1/|x - x'| e^(-i n t) dt
%
% it is assumed that x = (x,0,z) otherwise the integral above should pick up a
% phase factor out front of exp(i*n*phi), where phi is the azimuthal coordinate
% of x in cylindrical coordinates.
%
% The extra factor of rp (and maybe pi?) out front makes subsequent interfacing
% with RCIP slightly easier. Modes 0 through maxm are returned, with gval(1) =
% mode 0 and gval(maxm+1) = mode maxm. The function is even, so g_{-n} = g_n.
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
    xminus = (dr*dr + dz*dz)/2/r/r0;

    dxdr = (r^2 - r0^2 - (dz)^2)/2/r0/r^2;
    dxdz = 2*(dz)/2/r/r0;
    dxdr0 = (r0^2 - r^2 - (dz)^2)/2/r/r0^2;
    dxdz0 = -dxdz;


  
    %!
    %! if xminus is very small, use the forward recurrence
    %!

    iffwd = 0;
    if (x < 1.005)
        iffwd = 1;
        if ((x >= 1.0005d0) && (maxm > 163))
            iffwd = 0;
        end
        
        if ((x >= 1.00005d0) && (maxm > 503))
            iffwd = 0;
        end
        
        if ((x >= 1.000005d0) && (maxm > 1438))
            iffwd = 0;
        end
        
        if ((x >= 1.0000005d0) && (maxm > 4380))
            iffwd = 0;
        end
        
        if ((x >= 1.00000005d0) && (maxm > 12307))
            iffwd = 0;
        end
        
    end
    
    

    gvals = zeros(maxm+1,1);
    gdzs = zeros(maxm+1,1);
    gdrs = zeros(maxm+1,1);
    gdrps = zeros(maxm+1,1);
    

    if (iffwd == 1)

        [q0, q1, dq0] = chnk.axissymlap2d.qleg_half(xminus);
        dq1 = (-q0 + x*q1)/2/(x+1)/xminus;
        
        half = done/2;
    
        %fac = done/sqrt(r*r0)/4/pi^2;
        fac = 2*pi*sqrt(rp/r);
        gvals(1) = fac*q0;
        gvals(2) = fac*q1;

        derprev = fac*dq0;
        der = fac*dq1;
        
        % the z derivatives
        gdzs(1) = derprev*dxdz;
        gdzs(2) = der*dxdz;
        
        % the r derivatives
        gdrs(1) = fac*(dq0*dxdr - q0/2/r);
        gdrs(2) = fac*(dq1*dxdr - q1/2/r);

        % the rp derivatives
        gdrps(1) = fac*(dq0*dxdr0 - q0/2/r0);
        gdrps(2) = fac*(dq1*dxdr0 - q1/2/r0);

        % don't compute the zp derivatives, just minus the z derivatives
        %grad0s(2,0) = dq0*dxdz0
        %grad0s(2,1) = dq1*dxdz0

        % run upward recursion for the Q's to calculate them things
        for i = 1:(maxm-1)
            gvals(i+2) = (2*i*x*gvals(i+1) - (i-half)*gvals(i))/(i+half);
            dernext = (2*i*(gvals(i+1)+x*der) - (i-half)*derprev)/(i+half);
            gdrs(i+2) = (dernext*dxdr - gvals(i+2)/2/r);
            gdzs(i+2) = dernext*dxdz;
            gdrps(i+2) = (dernext*dxdr0 - gvals(i+2)/2/r0);
            derprev = der;
            der = dernext;
        end

        return
    end



    
    %!
    %! if here, xminus > .005, so run forward and backward recurrence
    %!

    %!
    %! run the recurrence, starting from maxm, until it has exploded
    %! for BOTH the values and derivatives
    %!
    done = 1;
    half = done/2;
    f = 1;
    fprev = 0;
    der = 1;
    derprev = 0;
    maxiter = 100000;
    upbound = 1.0e19;
  
    for i = maxm:maxiter
        fnext = (2*i*x*f - (i-half)*fprev)/(i+half);
        dernext = (2*i*(x*der+f) - (i-half)*derprev)/(i+half);
        if (abs(fnext) >= upbound) 
            if (abs(dernext) >= upbound)
                nterms = i+1;
                break
            end
        end
        fprev = f;
        f = fnext;
        derprev = der;
        der = dernext;
    end

    
    %!
    %! now start at nterms and recurse down to maxm
    %!

    if (nterms < 10)
        nterms = 10;
    end

    fnext = 0;
    f = 1;
    dernext = 0;
    der = 1;

    % run the downward recurrence
    for j = 1:(nterms-maxm+1)
        i = nterms-j+1;
        fprev = (2*i*x*f - (i+half)*fnext)/(i-half);
        fnext = f;
        f = fprev;
        derprev = (2*i*(x*der+f) - (i+half)*dernext)/(i-half);
        dernext = der;
        der = derprev;
    end

    gvals(maxm) = f;
    gvals(maxm+1) = fnext;

    ders = zeros(maxm+1,1);
    ders(maxm) = der;
    ders(maxm+1) = dernext;
  
    for j = 1:(maxm-1)
        i = maxm-1-j+1;
        gvals(i) = (2*i*x*gvals(i+1) - (i+half)*gvals(i+2))/(i-half);
        ders(i) = (2*i*(x*ders(i+1)+gvals(i+1)) - (i+half)*ders(i+2))/(i-half);
   end


    % !
    % ! normalize the values, and use a formula for the derivatives
    % !
    [q0, q1, dq0] = chnk.axissymlap2d.qleg_half(xminus);
    dq1 = (-q0 + x*q1)/2/(x+1)/xminus;
    % call axi_q2lege01(x, xminus, q0, q1, dq0, dq1)

    ratio = q0/gvals(1)*2*pi*sqrt(rp/r);
    
    for i = 1:(maxm+1)
        gvals(i) = gvals(i)*ratio;
    end

    ders(1) = dq0*2*pi*sqrt(rp/r);
    ders(2) = dq1*2*pi*sqrt(rp/r);
    for i = 2:maxm
       ders(i+1) = -(i-.5d0)*(gvals(i) - x*gvals(i+1))/(1+x)/xminus;
    end

    %
    % and scale the gradients properly everyone...
    %
    
    %fac = 1/sqrt(r*r0)/4/pi^2

    gdzs = ders*dxdz;
    gdrs = ders*dxdr - gvals/2/r;
    gdrps = ders*dxdr0 - gvals/2/r0;

    
    %for i = 1:(maxm+1)
    %   grads(1,i) = (ders(i)*dxdr - vals(i)/2/r)*fac
    % grads(2,i) = ders(i)*dxdz*fac
    %gdrs(i) = ders(i)*dxdr - gvals(i)/2/r;
        %gdzs(i) = ders(i)*dxdz;
    %   grad0s(1,i) = (ders(i)*dxdr0 - vals(i)/2/r0)*fac
    %   grad0s(2,i) = ders(i)*dxdz0*fac
    %    gvals(i) = gvals(i)*fac;
    %   vals(-i) = vals(i)
    %   grads(1,-i) = grads(1,i)
    %   grads(2,-i) = grads(2,i)
    %   grad0s(1,-i) = grad0s(1,i)
    %   grad0s(2,-i) = grad0s(2,i)
    %end


end
