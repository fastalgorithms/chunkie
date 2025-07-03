function submat = kern(srcinfo,targinfo,type,zk, vfun)
%CHNK.HELM2D_SHAPE_DER.KERN shape derivatives of 
% Laplace layer potential kernels in 2D
%
% Syntax: submat = kern(srcinfo, targinfo, type)

src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);

ivx = vfun(1);
ivy = vfun(2);
idvx = vfun(3);
idvy = vfun(4);
if strcmpi(type,'s')
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
    [val,grad] = chnk.helm2d.green(zk,src,targ);
    submat = vx_d.*grad(:,:,1)+vy_d.*grad(:,:,2);
    submat = submat + val.*(dvsx.*dx+dvsy.*dy)./(dx.*dx+dy.*dy);
end

if strcmpi(type,'sp')
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
    [~,grad,hess] = chnk.helm2d.green(zk,src,targ);
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
end

