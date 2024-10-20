function [gvals, gdzs, gdrs, gdrps] = g0funcall(r, rp, dr, z, zp, dz, maxm)
%
% chnk.axissymlap2d.g0funcall evaluates a collection of axisymmetric Laplace
% Green's functions, defined by the expression:
%
%     gfunc(n) = pi*rp * \int_0^{2\pi} 1/|x - x'| e^(-i n t) dt 
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


    %r = targ(1)
    %z = targ(2)
    r0 = rp;
    z0 = zp;
    rzero = sqrt(r*r + r0*r0 + dz*dz);
    alpha = 2*r*r0/rzero^2;
    x = 1/alpha;
    xminus = (dr*dr + dz*dz)/2/r/r0;

    dxdr = (r^2 - r0^2 - (dz)^2)/2/r0/r^2
    dxdz = 2*(dz)/2/r/r0
    dxdr0 = (r0^2 - r^2 - (dz)^2)/2/r/r0^2
    dxdz0 = -dxdz
  
    %!
    %! if xminus is very small, use the forward recurrence
    %!

    %!!!!print *, 'inside g0mall, x = ', x
  
  
    iffwd = 0
    if (x < 1.005)
        iffwd = 1
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

  
    if (iffwd = 1)

        %!!xminus = .000005d0
        %!!x = 1 + xminus
    
        call axi_q2lege01(x, xminus, q0, q1, dq0, dq1)

        done = 1
        half = done/2
        %pi = 4*atan(done)
    
        fac = done/sqrt(r*r0)/4/pi^2
        vals(0) = q0
        vals(1) = q1

        derprev = dq0
        der = dq1

        grads(1,0) = (dq0*dxdr - q0/2/r)
        grads(1,1) = (dq1*dxdr - q1/2/r)
        
        grads(2,0) = dq0*dxdz
        grads(2,1) = dq1*dxdz
        
        grad0s(1,0) = (dq0*dxdr0 - q0/2/r0)
        grad0s(1,1) = (dq1*dxdr0 - q1/2/r0)
        
        grad0s(2,0) = dq0*dxdz0
        grad0s(2,1) = dq1*dxdz0
        
        for i = 1:(maxm-1)
            vals(i+1) = (2*i*x*vals(i) - (i-half)*vals(i-1))/(i+half)
            dernext = (2*i*(vals(i)+x*der) - (i-half)*derprev)/(i+half)
            grads(1,i+1) = (dernext*dxdr - vals(i+1)/2/r)
            grads(2,i+1) = dernext*dxdz
            grad0s(1,i+1) = (dernext*dxdr0 - vals(i+1)/2/r0)
            grad0s(2,i+1) = dernext*dxdz0
            derprev = der
            der = dernext
            %! if (abs(vals(i+1)) > abs(vals(i))) then
            %!   print *
            %!   print *
            %!   call prin2('x = *', x, 1)
            %!   call prinf('i = *', i, 1)
            %!   call prin2('vals(i+1) = *', vals(i+1), 1)
            %!   call prin2('vals(i) = *', vals(i), 1)
            %!   stop
            %! endif
        end

        for i = 0:maxm
            vals(i) = vals(i)*fac
            grads(1,i) = grads(1,i)*fac
            grads(2,i) = grads(2,i)*fac
            grad0s(1,i) = grad0s(1,i)*fac
            grad0s(2,i) = grad0s(2,i)*fac
            vals(-i) = vals(i)
            grads(1,-i) = grads(1,i)
            grads(2,-i) = grads(2,i)
            grad0s(1,-i) = grad0s(1,i)
            grad0s(2,-i) = grad0s(2,i)
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
    done = 1
    half = done/2
    f = 1
    fprev = 0
    der = 1
    derprev = 0
    maxiter = 100000
    upbound = 1.0d19
  
    for i = maxm:maxiter
        fnext = (2*i*x*f - (i-half)*fprev)/(i+half)
        dernext = (2*i*(x*der+f) - (i-half)*derprev)/(i+half)
        if (abs(fnext) .ge. upbound) 
            if (abs(dernext) .ge. upbound)
                nterms = i+1
                exit
            end
        end
        fprev = f
        f = fnext
        derprev = der
        der = dernext
    end
    
    %!
    %! now start at nterms and recurse down to maxm
    %!

    if (nterms .lt. 10) then
        nterms = 10
    end

    fnext = 0
    f = 1
    dernext = 0
    der = 1


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
    % Make correct downwrd recurrence
    
  for i = nterm,maxm,-1
    fprev = (2*i*x*f - (i+half)*fnext)/(i-half)
    fnext = f
    f = fprev
    derprev = (2*i*(x*der+f) - (i+half)*dernext)/(i-half)
    dernext = der
    der = derprev
  enddo

  vals(maxm-1) = f
  vals(maxm) = fnext

  ders(maxm-1) = der
  ders(maxm) = dernext
  
  do i = maxm-1,1,-1
    vals(i-1) = (2*i*x*vals(i) - (i+half)*vals(i+1))/(i-half)
    ders(i-1) = (2*i*(x*ders(i)+vals(i)) &
        - (i+half)*ders(i+1))/(i-half)
  end do


  !
  ! normalize the values, and use a formula for the derivatives
  !
  call axi_q2lege01(x, xminus, q0, q1, dq0, dq1)

  done = 1
  pi = 4*atan(done)
  ratio = q0/vals(0)

  do i = 0,maxm
    vals(i) = vals(i)*ratio
  enddo

  ders(0) = dq0
  ders(1) = dq1
  do i = 2,maxm
    ders(i) = -(i-.5d0)*(vals(i-1) - x*vals(i))/(1+x)/xminus
  end do

  !
  ! and scale everyone...
  !
  fac = 1/sqrt(r*r0)/4/pi^2
  
  do i = 0,maxm
    grads(1,i) = (ders(i)*dxdr - vals(i)/2/r)*fac
    grads(2,i) = ders(i)*dxdz*fac
    grad0s(1,i) = (ders(i)*dxdr0 - vals(i)/2/r0)*fac
    grad0s(2,i) = ders(i)*dxdz0*fac
    vals(i) = vals(i)*fac
    vals(-i) = vals(i)
    grads(1,-i) = grads(1,i)
    grads(2,-i) = grads(2,i)
    grad0s(1,-i) = grad0s(1,i)
    grad0s(2,-i) = grad0s(2,i)
  end do
    




end
