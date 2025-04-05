function submat = kern(srcinfo,targinfo,type,zk,vfun)
%CHNK.HELM2DV.KERN standard Helmholtz layer potential kernels in 2D
%
% Syntax: submat = kern(srcinfo,targinfo,type,varargin)

src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);

if (nargin >4)
    if (isa(vfun,'function_handle'))
        if (strcmpi(type,'spdk') || strcmpi(type,'skv') || ...
                strcmpi(type,'sprimekv') || strcmpi(type,'dv') || ...
                strcmpi(type,'vSpp-Sppv') || ...
                strcmpi(type,'vSpp-Sppv+vdtauStau') || ...
                strcmpi(type,'dpv+sppv'))
            tpars = srcinfo.data(1,:);
            vvals = vfun(tpars);
        end
        if (strcmpi(type,'vsprime') || strcmpi(type,'vSpp-Sppv') || ...
                strcmpi(type,'vdtauStau') || ...
                strcmpi(type,'vSpp-Sppv+vdtauStau'))
            tpars = targinfo.data(1,:);
            vvalt = vfun(tpars);
            dt = targinfo.data(2,:);
            [~,dvals] = vfun(tpars);
            dvdt = dvals./dt;
        end
    elseif (isa(vfun,'double') && (numel(vfun)==2))
        ivv = vfun(1);
        ivd = vfun(2);
        if (strcmpi(type,'spdk') || strcmpi(type,'skv') || ...
                strcmpi(type,'sprimekv') || strcmpi(type,'dv') || ...
                strcmpi(type,'vSpp-Sppv') || ...
                strcmpi(type,'vSpp-Sppv+vdtauStau') || ...
                strcmpi(type,'dpv+sppv'))
            vvals = srcinfo.data(ivv,:);
        end
        if (strcmpi(type,'vsprime') || strcmpi(type,'vSpp-Sppv') || ...
                strcmpi(type,'vdtauStau') || ...
                strcmpi(type,'vSpp-Sppv+vdtauStau'))
            vvalt = targinfo.data(ivv,:);
            dvdt  = targinfo.data(ivd,:);
        end    
    end
end
if strcmpi(type,'spdk')

    xt = repmat(targinfo.r(1,:).',1,ns);
    yt = repmat(targinfo.r(2,:).',1,ns);
    xs = repmat(srcinfo.r(1,:),nt,1);
    ys = repmat(srcinfo.r(2,:),nt,1);

    targnorm = chnk.normal2d(targinfo);
    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);
    
    prefac = zk*((xt-xs).*nx+(yt-ys).*ny);

    submat = -prefac.*chnk.helm2d.green(zk,src,targ);
end

if strcmpi(type,'skv')
%    tpars = srcinfo.data(1,:);
%    vvals = vfun(tpars);
    kapps = chnk.curvature2d(srcinfo);
    vv = repmat(vvals,nt,1);
    ka = repmat(kapps,nt,1);
    submat = vv.*ka.*chnk.helm2d.green(zk,src,targ);
end

if strcmpi(type,'sprimekv')
%    tpars = srcinfo.data(1,:);
%    vvals = vfun(tpars);
    kapps = chnk.curvature2d(srcinfo);
    vv = repmat(vvals,nt,1);
    ka = repmat(kapps,nt,1);
    targnorm = chnk.normal2d(targinfo);
    [~,grad] = chnk.helm2d.green(zk,src,targ);
    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);
    submat = vv.*ka.*(grad(:,:,1).*nx + grad(:,:,2).*ny);
end

if strcmpi(type,'dv')
    %srcnorm = chnk.normal2d(srcinfo);
%    tpars = srcinfo.data(1,:);
%    vvals = vfun(tpars);
    [~,grad] = chnk.helm2d.green(zk,src,targ);
    vv = repmat(vvals,nt,1);
    submat = -vv.*(grad(:,:,1).*srcinfo.n(1,:) + grad(:,:,2).*srcinfo.n(2,:));
end

if strcmpi(type,'vsprime')

    targnorm = chnk.normal2d(targinfo);
    [~,grad] = chnk.helm2d.green(zk,src,targ);
    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);
    vv = repmat(vvalt.',1,ns);
    submat = vv.*(grad(:,:,1).*nx + grad(:,:,2).*ny);
end

if strcmpi(type,'vSpp-Sppv')
    targnorm = targinfo.n;
    srcnorm = srcinfo.n;
    [~,~,hess] = chnk.helm2d.green(zk,src,targ);



    vt = repmat(vvalt.',1,ns);
    vs = repmat(vvals,nt,1);

    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
    submat = (hess(:,:,1).*nxtarg.*nxtarg + hess(:,:,2).*(nytarg.*nxtarg+nxtarg.*nytarg)...
        + hess(:,:,3).*nytarg.*nytarg);
    submat = submat.*(vt-vs);
end

if strcmpi(type,'vdtauStau')

%    tpars = targinfo.data(1,:);
%    dt = targinfo.data(2,:);
%    [~,dvals] = vfun(tpars);
%    dvdt = dvals./dt;
    %dvdt = dvals;
    vv = repmat(dvdt.',1,ns);

    targnorm = chnk.normal2d(targinfo);
    [~,grad] = chnk.helm2d.green(zk,src,targ);
    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);

    submat = vv.*(-grad(:,:,1).*ny + grad(:,:,2).*nx);
end

if strcmpi(type,'vSpp-Sppv+vdtauStau')
    targnorm = targinfo.n;
    srcnorm = srcinfo.n;
    [~,~,hess] = chnk.helm2d.green(zk,src,targ);

%    tparst = targinfo.data(1,:);
%    vvalst = vfun(tparst);

%    tparss = srcinfo.data(1,:);
%    vvalss = vfun(tparss);

    vt = repmat(vvalt.',1,ns);
    vs = repmat(vvals,nt,1);

    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
    submat = (hess(:,:,1).*nxtarg.*nxtarg + hess(:,:,2).*(nytarg.*nxtarg+nxtarg.*nytarg)...
        + hess(:,:,3).*nytarg.*nytarg);
    submat = submat.*(vt-vs);

%    tpars = targinfo.data(1,:);
%    dt = targinfo.data(2,:);
%    [~,dvals] = vfun(tpars);
%    dvdt = dvals./dt;
    vv = repmat(dvdt.',1,ns);

    targnorm = chnk.normal2d(targinfo);
    [~,grad] = chnk.helm2d.green(zk,src,targ);
    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);

    submat2 = vv.*(-grad(:,:,1).*ny + grad(:,:,2).*nx);

    submat = submat - submat2;

end

if strcmpi(type,'dpv+sppv')
%    tpars = srcinfo.data(1,:);
%    vvals = vfun(tpars);
    vv = repmat(vvals,nt,1);
    targnorm = targinfo.n;
    srcnorm = srcinfo.n;
    [~,~,hess] = chnk.helm2d.green(zk,src,targ);
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



