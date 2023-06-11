
function submat= kern(zk,srcinfo,targinfo,type,varargin)
%CHNK.HELM2D.KERN standard Helmholtz layer potential kernels in 2D
% 
% Syntax: submat = chnk.heml2d.kern(zk,srcinfo,targingo,type,varargin)
%
% Let x be targets and y be sources for these formulas, with
% n_x and n_y the corresponding unit normals at those points
% (if defined). Note that the normal information is obtained
% by taking the perpendicular to the provided tangential deriviative
% info and normalizing  
%  
% Kernels based on G(x,y) = i/4 H_0^{(1)}(zk |x-y|)
%
% D(x,y) = \nabla_{n_y} G(x,y)
% S(x,y) = G(x,y)
% S'(x,y) = \nabla_{n_x} G(x,y)
%
% Input:
%   zk - complex number, Helmholtz wave number
%   srcinfo - description of sources in ptinfo struct format, i.e.
%                ptinfo.r - positions (2,:) array
%                ptinfo.d - first derivative in underlying
%                     parameterization (2,:)
%                ptinfo.d2 - second derivative in underlying
%                     parameterization (2,:)
%   targinfo - description of targets in ptinfo struct format,
%                if info not relevant (d/d2) it doesn't need to
%                be provided. sprime requires tangent info in
%                targinfo.d
%   type - string, determines kernel type
%                type == 'd', double layer kernel D
%                type == 's', single layer kernel S
%                type == 'sprime', normal derivative of single
%                      layer S'
%                type == 'c', combined layer kernel D + i eta S
%   varargin{1} - eta in the combined layer formula, otherwise
%                does nothing
%
% Output:
%   submat - the evaluation of the selected kernel for the
%            provided sources and targets. the number of
%            rows equals the number of targets and the
%            number of columns equals the number of sources  
%
% see also CHNK.HELM2D.GREEN
  
src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);

if strcmpi(type,'d')
  srcnorm = srcinfo.n;
  [~,grad] = chnk.helm2d.green(zk,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  submat = -(grad(:,:,1).*nx + grad(:,:,2).*ny);
end

if strcmpi(type,'sprime')
  targnorm = targinfo.n;
  [~,grad] = chnk.helm2d.green(zk,src,targ);
  nx = repmat((targnorm(1,:)).',1,ns);
  ny = repmat((targnorm(2,:)).',1,ns);

  submat = (grad(:,:,1).*nx + grad(:,:,2).*ny);
end

if strcmpi(type,'sdtau')
  targtan = targinfo.d;
  [~,grad] = chnk.helm2d.green(zk,src,targ);
  dx = repmat((targtan(1,:)).',1,ns);
  dy = repmat((targtan(2,:)).',1,ns);
  dn = sqrt(dx.*dx+dy.*dy);
  submat = (grad(:,:,1).*dx./dn + grad(:,:,2).*dy)./dn;
end

if strcmpi(type,'s')
  submat = chnk.helm2d.green(zk,src,targ);
end

if strcmpi(type,'dprime')
  targnorm = targinfo.n;
  srcnorm = srcinfo.n;
  [~,~,hess] = chnk.helm2d.green(zk,src,targ);
  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);
  submat = -(hess(:,:,1).*nxsrc.*nxtarg + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
end

if strcmpi(type,'c')
  srcnorm = srcinfo.n;
  coef = varargin{1};
  [submats,grad] = chnk.helm2d.green(zk,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  submatd = -(grad(:,:,1).*nx + grad(:,:,2).*ny);
  submat = coef(1)*submatd + coef(2)*submats;
end

if strcmpi(type,'all')
  %targnorm = chnk.normal2d(targinfo);
  %srcnorm = chnk.normal2d(srcinfo);
  % added by Shidong Jiang to avoid O(N^2) calculation of normals
  targnorm = targinfo.n;
  srcnorm = srcinfo.n;
  cc = varargin{1};
  
  submat = zeros(2*nt,2*ns);
  % S
  [submats,grad,hess] = chnk.helm2d.green(zk,src,targ);
  %[submat(1:2:2*nt,2:2:2*ns),grad,hess] = chnk.helm2d.green(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);
  % D'
  %submat(2:2:2*nt,1:2:2*ns) = -(hess(:,:,1).*nxsrc.*nxtarg ... 
  submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
      + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
  % S'
%  submat(2:2:2*nt,2:2:2*ns) = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
  submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
  % D
%  submat(1:2:2*nt,1:2:2*ns)  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);
  submatd  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);
  
  
  submat(1:2:2*nt,1:2:2*ns) = submatd*cc(1,1);
  submat(1:2:2*nt,2:2:2*ns) = submats*cc(1,2);
  submat(2:2:2*nt,1:2:2*ns) = submatdp*cc(2,1);
  submat(2:2:2*nt,2:2:2*ns) = submatsp*cc(2,2);
end

if strcmpi(type,'eval trans')

  srcnorm = srcinfo.n;
  [submats,grad] = chnk.helm2d.green(zk,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  submatd = -(grad(:,:,1).*nx + grad(:,:,2).*ny);

  submat = zeros(1*nt,2*ns);
  submat(1:1:1*nt,1:2:2*ns) = submatd;
  submat(1:1:1*nt,2:2:2*ns) = submats;
end

if strcmpi(type,'c and cprime')
  targnorm = targinfo.n;
  srcnorm = srcinfo.n;
  coef = varargin{1};

  submat = zeros(2*nt,1*ns);
  % S
  [submats,grad,hess] = chnk.helm2d.green(zk,src,targ);
  %[submat(1:2:2*nt,2:2:2*ns),grad,hess] = chnk.helm2d.green(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);
  % D'
  %submat(2:2:2*nt,1:2:2*ns) = -(hess(:,:,1).*nxsrc.*nxtarg ... 
  submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
      + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
  % S'
%  submat(2:2:2*nt,2:2:2*ns) = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
  submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
  % D
%  submat(1:2:2*nt,1:2:2*ns)  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);
  submatd  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);
  
  submatc      = coef(1,1) * submatd  + coef(1,2) * submats;
  submatcprime = coef(2,1) * submatdp + coef(2,2) * submatsp;
  
  submat(1:2:2*nt,1:1:1*ns) = submatc;
  submat(2:2:2*nt,1:1:1*ns) = submatcprime;
end

if strcmpi(type,'eval')
  coef = varargin{1};
  
  srcnorm = srcinfo.n;
  
  submat = zeros(nt,2*ns);
  % S
  [submats,grad] = chnk.helm2d.green(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  % D
  submatd  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);
    
  submat(:,1:2:2*ns) = coef(1)*submatd;
  submat(:,2:2:2*ns) = coef(2)*submats;
end


if strcmpi(type,'evalg')
  coef = varargin{1};
  
  srcnorm = srcinfo.n;
  
  submat = zeros(nt,ns,6);
  % S
  [submats,grad,hess] = chnk.helm2d.green(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  % D
  submatd  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);
    
  submat(:,:,1) = coef*submatd;
  submat(:,:,2) = submats;
  submat(:,:,3) = -coef*(hess(:,:,1).*nxsrc + hess(:,:,2).*nysrc);
  submat(:,:,4) = grad(:,:,1);
  submat(:,:,5) = -coef*(hess(:,:,2).*nxsrc + hess(:,:,3).*nysrc);
  submat(:,:,6) = grad(:,:,2);
end



if strcmpi(type,'trans1')
  %targnorm = chnk.normal2d(targinfo);
  %srcnorm = chnk.normal2d(srcinfo);
  % added by Shidong Jiang to avoid O(N^2) calculation of normals
  targnorm = targinfo.n;
  srcnorm = srcinfo.n;
  
  
  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);

  
  submat = zeros(2*nt,2*ns);

  xs = repmat(src(1,:),nt,1);
  ys = repmat(src(2,:),nt,1);

  [m,n] = size(xs);
  xt = repmat(targ(1,:).',1,ns);
  yt = repmat(targ(2,:).',1,ns);

  rx = xt-xs;
  ry = yt-ys;

  rx2 = rx.*rx;
  ry2 = ry.*ry;

  r2 = rx2+ry2;

  r = sqrt(r2);
  k1 = repmat(targinfo.data(1,:).',1,ns);
  k2 = repmat(targinfo.data(2,:).',1,ns);
  c1 = repmat(targinfo.data(3,:).',1,ns);
  c2 = repmat(targinfo.data(4,:).',1,ns);
  
  harg1  = k1.*r;
  hargsz = size(harg1);
  harg2  = k2.*r;
  hargsz2= size(harg2);
  harg = [harg1(:),harg2(:)];
  %szharg = size(harg)
  [h00,h11] = chnk.helm2d.besselh01(harg);
   %[h00,h11] = hankm103(harg);
  %szh00 = size(h00)
  h0 = reshape(h00(:,1),hargsz);
  h1 = reshape(h11(:,1),hargsz);
  %[h0,h1] = hankm103(k1.*r);
  %h0 = besselh(0,1,k1.*r);
  submats1 = 0.25*1i*h0;
  
  h1 = h1./r;
  %h1 = besselh(1,1,k1.*r)./r;
  
  grad1 = zeros(m,n,2);
    
  ck = 0.25*1i*k1;
  grad1(:,:,1) = -ck.*h1.*rx;
  grad1(:,:,2) = -ck.*h1.*ry;    

  hess1 = zeros(m,n,3);

  h2 = 2*h1./k1-h0;
  tmp1 = (rx2-ry2).*h1./r2;
  tmp2 = k1.*h0./r2;
    
    
  hess1(:,:,1) = ck.*(tmp1 - rx2.*tmp2);
  hess1(:,:,2) = ck.*k1.*rx.*ry.*h2./r2;
    
  hess1(:,:,3) = ck.*(-tmp1 - ry2.*tmp2);

  %[submat(1:2:2*nt,2:2:2*ns),grad,hess] = chnk.helm2d.green(zk,src,targ);
  % D'
  %submat(2:2:2*nt,1:2:2*ns) = -(hess(:,:,1).*nxsrc.*nxtarg ... 
  submatdp1 = -(hess1(:,:,1).*nxsrc.*nxtarg ...
      + hess1(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess1(:,:,3).*nysrc.*nytarg);
  % S'
%  submat(2:2:2*nt,2:2:2*ns) = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
  submatsp1 = (grad1(:,:,1).*nxtarg + grad1(:,:,2).*nytarg);
  % D
%  submat(1:2:2*nt,1:2:2*ns)  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);
  submatd1  = -(grad1(:,:,1).*nxsrc + grad1(:,:,2).*nysrc);
 
  h0 = reshape(h00(:,2),hargsz2);
  h1 = reshape(h11(:,2),hargsz2);
  %[h0,h1] = hankm103(k2.*r);
  %h0 = besselh(0,1,k2.*r);
  submats2 = 0.25*1i*h0;
  
  %h1 = besselh(1,1,k2.*r)./r;
  h1 = h1./r;
  
  grad2 = zeros(m,n,2);
    
  ck = 0.25*1i*k2;
  grad2(:,:,1) = -ck.*h1.*rx;
  grad2(:,:,2) = -ck.*h1.*ry;    

  hess2 = zeros(m,n,3);

  h2 = 2*h1./k2-h0;
  tmp1 = (rx2-ry2).*h1./r2;
  tmp2 = k2.*h0./r2;
    
    
  hess2(:,:,1) = ck.*(tmp1 - rx2.*tmp2);
  hess2(:,:,2) = ck.*k2.*rx.*ry.*h2./r2;
    
  hess2(:,:,3) = ck.*(-tmp1 - ry2.*tmp2);

  %[submat(1:2:2*nt,2:2:2*ns),grad,hess] = chnk.helm2d.green(zk,src,targ);
  % D'
  %submat(2:2:2*nt,1:2:2*ns) = -(hess(:,:,1).*nxsrc.*nxtarg ... 
  submatdp2 = -(hess2(:,:,1).*nxsrc.*nxtarg ...
      + hess2(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess2(:,:,3).*nysrc.*nytarg);
  % S'
%  submat(2:2:2*nt,2:2:2*ns) = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
  submatsp2 = (grad2(:,:,1).*nxtarg + grad2(:,:,2).*nytarg);
  % D
%  submat(1:2:2*nt,1:2:2*ns)  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);
  submatd2  = -(grad2(:,:,1).*nxsrc + grad2(:,:,2).*nysrc);
  
  
  alpha1 = 2./(c1+c2);
  alpha2 = 2./(1./c1+1./c2);
    
  
  submat(1:2:2*nt,1:2:2*ns) = alpha1.*(c2.*submatd2-c1.*submatd1);
  submat(1:2:2*nt,2:2:2*ns) = alpha1.*(submats2-submats1);
  submat(2:2:2*nt,1:2:2*ns) = -alpha2.*(submatdp2-submatdp1);
  submat(2:2:2*nt,2:2:2*ns) = -alpha2.*(1./c2.*submatsp2-1./c1.*submatsp1);
end
