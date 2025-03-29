function submat = kern(srcinfo,targinfo,type,vfun)
%CHNK.LAP2D.KERN standard Laplace layer potential kernels in 2D
%
% Syntax: submat = kern(srcinfo,targinfo,type,varargin)

src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);

if strcmpi(type,'skv')
    tpars = srcinfo.data(1,:);
    vvals = vfun(tpars);
    kapps = chnk.curvature2d(srcinfo);
    vv = repmat(vvals,nt,1);
    ka = repmat(kapps,nt,1);
    submat = vv.*ka.*chnk.lap2d.green(src,targ);
end

if strcmpi(type,'sprimekv')
    tpars = srcinfo.data(1,:);
    vvals = vfun(tpars);
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
    tpars = srcinfo.data(1,:);
    vvals = vfun(tpars);
    [~,grad] = chnk.lap2d.green(src,targ,true);
    vv = repmat(vvals,nt,1);
    submat = -vv.*(grad(:,:,1).*srcinfo.n(1,:) + grad(:,:,2).*srcinfo.n(2,:));
end

if strcmpi(type,'vsprime')

    tpars = targinfo.data(1,:);
    vvals = vfun(tpars);
    targnorm = chnk.normal2d(targinfo);
    [~,grad] = chnk.lap2d.green(src,targ,true);
    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);
    vv = repmat(vvals.',1,ns);
    submat = vv.*(grad(:,:,1).*nx + grad(:,:,2).*ny);
end

if strcmpi(type,'vSpp-Sppv')
    targnorm = targinfo.n;
    srcnorm = srcinfo.n;
    [~,~,hess] = chnk.lap2d.green(src,targ);

    tparst = targinfo.data(1,:);
    vvalst = vfun(tparst);

    tparss = srcinfo.data(1,:);
    vvalss = vfun(tparss);

    vt = repmat(vvalst.',1,ns);
    vs = repmat(vvalss,nt,1);

    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
    submat = (hess(:,:,1).*nxtarg.*nxtarg + hess(:,:,2).*(nytarg.*nxtarg+nxtarg.*nytarg)...
        + hess(:,:,3).*nytarg.*nytarg);
    submat = submat.*(vt-vs);
end

if strcmpi(type,'vdtauStau')

    tpars = targinfo.data(1,:);
    dt = targinfo.data(2,:);
    [~,dvals] = vfun(tpars);
    dvdt = dvals./dt;
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

    tparst = targinfo.data(1,:);
    vvalst = vfun(tparst);

    tparss = srcinfo.data(1,:);
    vvalss = vfun(tparss);

    vt = repmat(vvalst.',1,ns);
    vs = repmat(vvalss,nt,1);

    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
    submat = (hess(:,:,1).*nxtarg.*nxtarg + hess(:,:,2).*(nytarg.*nxtarg+nxtarg.*nytarg)...
        + hess(:,:,3).*nytarg.*nytarg);
    submat = submat.*(vt-vs);

    tpars = targinfo.data(1,:);
    dt = targinfo.data(2,:);
    [~,dvals] = vfun(tpars);
    dvdt = dvals./dt;
    vv = repmat(dvdt.',1,ns);

    targnorm = chnk.normal2d(targinfo);
    [~,grad] = chnk.lap2d.green(src,targ,true);
    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);

    submat2 = vv.*(-grad(:,:,1).*ny + grad(:,:,2).*nx);

    submat = submat - submat2;

end

if strcmpi(type,'dpv+sppv')
    tpars = srcinfo.data(1,:);
    vvals = vfun(tpars);
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


