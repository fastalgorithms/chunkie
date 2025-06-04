function submat = kern(srcinfo,targinfo,type,vfun)
%CHNK.LAP2D.KERN standard Laplace layer potential kernels in 2D
%
% Syntax: submat = kern(srcinfo,targinfo,type,varargin)

src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);

ivx = vfun(1);
ivy = vfun(2);
idvx = vfun(3);
idvy = vfun(4);
if strcmpi(type,'dsdz')
    vtx = repmat(targinfo.data(ivx,:).',1,ns);
    vty = repmat(targinfo.data(ivy,:).',1,ns);
    vsx = repmat(srcinfo.data(ivx,:)   ,nt,1);
    vsy = repmat(srcinfo.data(ivy,:)   ,nt,1);
    dvsx = repmat(srcinfo.data(idvx,:)   ,nt,1);
    dvsy = repmat(srcinfo.data(idvy,:)   ,nt,1);
    dx = repmat(srcinfo.d(1,:),nt,1);
    dy = repmat(srcinfo.d(2,:),nt,1);
    vx_d = vtx-vsx;
    vy_d = vty-vsy;
    [val,grad] = chnk.lap2d.green(src,targ);
    submat = vx_d.*grad(:,:,1)+vy_d.*grad(:,:,2);
    submat = submat + val.*(dvsx.*dx+dvsy.*dy)./(dx.*dx+dy.*dy);
end

if strcmpi(type,'dspdz')
    ntx = repmat(targinfo.n(1,:).',1,ns);
    nty = repmat(targinfo.n(2,:).',1,ns);
    vtx = repmat(targinfo.data(ivx,:).',1,ns);
    vty = repmat(targinfo.data(ivy,:).',1,ns);
    vsx = repmat(srcinfo.data(ivx,:)   ,nt,1);
    vsy = repmat(srcinfo.data(ivy,:)   ,nt,1);
    dvsx = repmat(srcinfo.data(idvx,:)   ,nt,1);
    dvsy = repmat(srcinfo.data(idvy,:)   ,nt,1);
    dvtx = repmat(targinfo.data(idvx,:).',1,ns);
    dvty = repmat(targinfo.data(idvy,:).',1,ns);
    dx = repmat(srcinfo.d(1,:),nt,1);
    dy = repmat(srcinfo.d(2,:),nt,1);
    dxt = repmat(targinfo.d(1,:).',1,ns);
    dyt = repmat(targinfo.d(2,:).',1,ns);
    vx_d = vtx-vsx;
    vy_d = vty-vsy;
    [~,grad,hess] = chnk.lap2d.green(src,targ);
    kp = grad(:,:,1).*ntx+grad(:,:,2).*nty;

    a11 = ntx.*vx_d;
    a12 = ntx.*vy_d;
    a21 = nty.*vx_d;
    a22 = nty.*vy_d;

    ds  = sqrt(dxt.^2+dyt.^2);
    prefac = (dvtx.*dxt+dvty.*dyt)./ds.^2;
    b1 = dvty./ds - ntx.*prefac;
    b2 =-dvtx./ds - nty.*prefac;

    t1 = hess(:,:,1).*a11+hess(:,:,2).*a12+hess(:,:,2).*a21+...
        hess(:,:,3).*a22;
    t2 = b1.*grad(:,:,1)+b2.*grad(:,:,2);
    t3 = kp.*(dvsx.*dx+dvsy.*dy)./(dx.*dx+dy.*dy);

    submat = t1 + t2 + t3;
    %submat = vx_d.*grad(:,:,1)+vy_d.*grad(:,:,2);
    %submat = submat + val.*(dvsx.*dx+dvsy.*dy)./(dx.*dx+dy.*dy);
end

if strcmpi(type,'skv')
%    tpars = srcinfo.data(1,:);
%    vvals = vfun(tpars);
    kapps = chnk.curvature2d(srcinfo);
    vv = repmat(vvals,nt,1);
    ka = repmat(kapps,nt,1);
    submat = vv.*ka.*chnk.lap2d.green(src,targ);
end

if strcmpi(type,'sprimekv')
%    tpars = srcinfo.data(1,:);
%    vvals = vfun(tpars);
    kapps = chnk.curvature2d(srcinfo);
    vv = repmat(vvals,nt,1);
    ka = repmat(kapps,nt,1);
    targnorm = chnk.normal2d(targinfo);
    [~,grad] = chnk.lap2d.green(src,targ,true);
    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);
    submat = vv.*ka.*(grad(:,:,1).*nx + grad(:,:,2).*ny);
end

if strcmpi(type,'dv')
    %srcnorm = chnk.normal2d(srcinfo);
    %tpars = srcinfo.data(1,:);
    %vvals = vfun(tpars);
    [~,grad] = chnk.lap2d.green(src,targ,true);
    vv = repmat(vvals,nt,1);
    submat = -vv.*(grad(:,:,1).*srcinfo.n(1,:) + grad(:,:,2).*srcinfo.n(2,:));
end

if strcmpi(type,'vsprime')

    %tpars = targinfo.data(1,:);
    %vvals = vfun(tpars);
    targnorm = chnk.normal2d(targinfo);
    [~,grad] = chnk.lap2d.green(src,targ,true);
    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);
    vv = repmat(vvalt.',1,ns);
    submat = vv.*(grad(:,:,1).*nx + grad(:,:,2).*ny);
end

if strcmpi(type,'vSpp-Sppv')
    targnorm = targinfo.n;
    srcnorm = srcinfo.n;
    [~,~,hess] = chnk.lap2d.green(src,targ);

    %tparst = targinfo.data(1,:);
    %vvalst = vfun(tparst);

    %tparss = srcinfo.data(1,:);
    %vvalss = vfun(tparss);

    vt = repmat(vvalt.',1,ns);
    vs = repmat(vvals,nt,1);

    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
    submat = (hess(:,:,1).*nxtarg.*nxtarg + hess(:,:,2).*(nytarg.*nxtarg+nxtarg.*nytarg)...
        + hess(:,:,3).*nytarg.*nytarg);
    submat = submat.*(vt-vs);
end

if strcmpi(type,'vdtauStau')

    %tpars = targinfo.data(1,:);
    %dt = targinfo.data(2,:);
    %[~,dvals] = vfun(tpars);
    %dvdt = dvals./dt;
    %dvdt = dvals;
    vv = repmat(dvdt.',1,ns);

    targnorm = chnk.normal2d(targinfo);
    [~,grad] = chnk.lap2d.green(src,targ,true);
    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);

    submat = vv.*(-grad(:,:,1).*ny + grad(:,:,2).*nx);
end

if strcmpi(type,'vSpp-Sppv+vdtauStau')
    targnorm = targinfo.n;
    srcnorm = srcinfo.n;
    [~,~,hess] = chnk.lap2d.green(src,targ);

    %tparst = targinfo.data(1,:);
    %vvalst = vfun(tparst);

    %tparss = srcinfo.data(1,:);
    %vvalss = vfun(tparss);

    vt = repmat(vvalt.',1,ns);
    vs = repmat(vvals,nt,1);

    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
    submat = (hess(:,:,1).*nxtarg.*nxtarg + hess(:,:,2).*(nytarg.*nxtarg+nxtarg.*nytarg)...
        + hess(:,:,3).*nytarg.*nytarg);
    submat = submat.*(vt-vs);

    %tpars = targinfo.data(1,:);
    %dt = targinfo.data(2,:);
    %[~,dvals] = vfun(tpars);
    %dvdt = dvals./dt;
    vv = repmat(dvdt.',1,ns);

    targnorm = chnk.normal2d(targinfo);
    [~,grad] = chnk.lap2d.green(src,targ,true);
    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);

    submat2 = vv.*(-grad(:,:,1).*ny + grad(:,:,2).*nx);

    submat = submat - submat2;

end

if strcmpi(type,'dpv+sppv')
    %tpars = srcinfo.data(1,:);
    %vvals = vfun(tpars);
    vv = repmat(vvals,nt,1);
    targnorm = targinfo.n;
    srcnorm = srcinfo.n;
    [~,~,hess] = chnk.lap2d.green(src,targ);
    nxsrc = repmat(srcnorm(1,:),nt,1);
    nysrc = repmat(srcnorm(2,:),nt,1);
    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
    submat = -(hess(:,:,1).*nxsrc.*nxtarg + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
        + hess(:,:,3).*nysrc.*nytarg);
    submat2= (hess(:,:,1).*nxtarg.*nxtarg + hess(:,:,2).*(nytarg.*nxtarg+nxtarg.*nytarg)...
        + hess(:,:,3).*nytarg.*nytarg);
    submat = (submat+submat2).*vv;
end

if strcmpi(type,'stau')
    targnorm = chnk.normal2d(targinfo);
    [~,grad] = chnk.lap2d.green(src,targ,true);
    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);

    submat = (-grad(:,:,1).*ny + grad(:,:,2).*nx);
end

if strcmpi(type,'sgrad')
    [~,grad] = chnk.lap2d.green(src,targ,true);
    submat = reshape(permute(grad,[3,1,2]),2*nt,ns);
end

if strcmpi(type,'dgrad')
    [~,~,hess] = chnk.lap2d.green(src,targ,true);
    submat = -(hess(:,:,1:2).*srcinfo.n(1,:)+hess(:,:,2:3).*srcinfo.n(2,:));
    submat = reshape(permute(submat,[3,1,2]),2*nt,ns);
end

if strcmpi(type,'dprime')
    targnorm = targinfo.n;
    srcnorm = srcinfo.n;
    [~,~,hess] = chnk.lap2d.green(src,targ);
    nxsrc = repmat(srcnorm(1,:),nt,1);
    nysrc = repmat(srcnorm(2,:),nt,1);
    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
    submat = -(hess(:,:,1).*nxsrc.*nxtarg + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
        + hess(:,:,3).*nysrc.*nytarg);
end


