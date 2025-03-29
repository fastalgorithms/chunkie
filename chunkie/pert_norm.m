function [r,d] = pert_norm(t,func,vfunc,s)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   this is just a convenient helper routine which enables
    %   the construction of normal perturbations of a given shape.
    %
    %   func - the given parameterization:
    %       [r,d,d2] = func(t)
    %   vfunc - a function handle giving the coefficient of the normal 
    %       in the parameterization:
    %       [v,dvdt] = vfunc(t)
    %   s - a scalar, multiplying the perturbation
    %
    %   In particular, r = func(t) + s*v(t)*n(t)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sz = size(t);
    t  = t(:).';
    [fr,fd,fd2] = func(t);
    
    fdn = sqrt(fd(1,:).^2+fd(2,:).^2);
    taux = fd(1,:)./fdn;
    tauy = fd(2,:)./fdn;

    fdnt = (fd(1,:).*fd2(1,:)+fd(2,:).*fd2(2,:))./fdn;
    dtaux = fd2(1,:)./fdn - fd(1,:)./(fdn.^2).*fdnt;
    dtauy = fd2(2,:)./fdn - fd(2,:)./(fdn.^2).*fdnt;

    nx   = tauy;
    ny   =-taux;
    dnx  = dtauy;
    dny  =-dtaux;


    [v,dv] = vfunc(t);
    r = fr;
    r(1,:) = r(1,:) + nx.*v*s;
    r(2,:) = r(2,:) + ny.*v*s;
    
    d = fd;
    d(1,:) = d(1,:) + nx.*dv*s + dnx.*v*s;
    d(2,:) = d(2,:) + ny.*dv*s + dny.*v*s;

    r = reshape(r,[2,squeeze(sz)]);
    d = reshape(d,[2,squeeze(sz)]);
end