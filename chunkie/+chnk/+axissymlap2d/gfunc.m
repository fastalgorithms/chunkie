function [gval,gdz,gdr,gdrp,gdzz,gdrrp,gdrz,gdrpz] = gfunc(r,rp,dr,dz,n)
    % chnk.axissymlap2d.gfunc
    % green's function kernel normalized by area
    domega = 1/(2*pi);
    %wn = @(n) 2*pi^(n/2)/gamma(n/2);
    t0 = (r == 0);
    s0 = (rp == 0);
    st0 = ~((~t0).*(~s0));
    dz2 = dz.^2;
    gval0 = 1/2*rp./(sqrt(r.^2+rp.^2+dz2));
    gdrp0 = 1/2*rp.^2./(sqrt(rp.^2+dz.^2).^3); 
    gdz0  = 1/2*rp.*dz./(sqrt(rp.^2+dz.^2).^3);

    t = (dz.^2+dr.^2)./(2.*r.*rp);
    factor = (rp./r).^((n-2)/2);
    %[q0,~,q0d,q0dd] = chnk.axissymlap2d.qleg_half(t); % TODO: q0dd for D'
    %[~,~,~,q0dd] = chnk.axissymlap2d.qleg_half(t); % TODO: q0dd for D'
    inear = (t<1e-2);
    xm1n = t(inear);
    xm1f = t(~inear);
    [q0f,~,qdf,qddf] = chnk.axissymlap2d.runbackward(xm1f,n); 
    [q0n,~,qdn,qddn] = chnk.axissymlap2d.runforward(xm1n,n); 
    q0 = zeros(size(t)); % q0 = Q_{n/2-2}
    q0d = q0; q0dd = q0;
    q0(inear) = q0n; q0(~inear) = q0f;
    q0d(inear) = qdn; q0d(~inear) = qdf;
    q0dd(inear) = qddn; q0dd(~inear) = qddf;

    gval = domega*factor.*q0;
    gdz  =-domega*factor.*q0d ...
           ./(rp.*r).*dz;
    rpfac = -r*(n-2)/2.*q0+(-(1+t).*r+rp).*q0d;
    gdrp  = -domega*factor./(rp.*r) ...
          .*rpfac;
    rfac = -rp*(n-2)/2.*q0+(-(1+t).*rp+r).*q0d;
    gdr  = -domega*factor./(rp.*r) ...
          .*rfac;

    gval(st0) = gval0(st0);
    gdr(st0) = 0;
    gdrp(st0) = gdrp0(st0);
    gdz(st0) = gdz0(st0);
     
    % t derivatives * rrp
    tdr = r-rp.*(1+t);
    tdrp = rp-(1+t).*r;
    tdz = -dz;
    rrp = rp.*r;

    % compute for hessian
    gdzz = domega*factor.*(q0dd.*(dz./rrp).^2+q0d./rrp); 

    %q0dterm = 2*(1+t)-(r./rp+rp./r)*3/2;
    q0dterm = (n-1)*(1+t)-(r./rp+rp./r)*n/2;
    gdrrp = (n-2)^2*q0/4 + q0d.*q0dterm + q0dd./rrp.*tdr.*tdrp;
    gdrrp = -gdrrp*domega.*factor./rrp;
    
    gdrz = q0dd.*tdr - q0d.*rp*n/2;
    gdrz = -gdrz*domega.*factor.*dz./rrp.^2; 

    gdrpz = q0dd.*tdrp - q0d.*r*n/2;
    gdrpz = -gdrpz*domega.*factor.*dz./rrp.^2; 
    
end