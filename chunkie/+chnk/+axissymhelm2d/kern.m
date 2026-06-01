function submat = kern(zks, srcinfo, targinfo, origin, type, varargin)
%CHNK.AXISSYMHELM2D.KERN axissymmetric Helmholtz layer potential kernels in 2D
% 
% Syntax: submat = chnk.axissymhelm2d.kern(zk,srcinfo,targingo,type,htables)
%
% Let x be targets and y be sources for these formulas, with
% n_x and n_y the corresponding unit normals at those points
% (if defined). Note that the normal information is obtained
% by taking the perpendicular to the provided tangential deriviative
% info and normalizing  
%
% Here the first and second components correspond to the r and z
% coordinates respectively. 
%
% Kernels based on G(x,y) = \int_{0}^{\pi} e^{i k d(t)}/(d(t)) \, dt \, 
% where d(t) = \sqrt(r^2 + r'^2 - 2rr' \cos(t) + (z-z')^2) with
% x = (r,z), and y = (r',z')
%
% D(x,y) = \nabla_{n_y} G(x,y)
% S(x,y) = G(x,y)
% S'(x,y) = \nabla_{n_x} G(x,y)
% D'(x,y) = \nabla_{n_x} \nabla_{n_y} G(x,y)
%
% Input:
%   zk - complex number, Helmholtz wave number, must be purely real, or purely 
%        imaginary, doesn't support other complex wavenumbers yet
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
%                type == 'dprime', normal derivative of double layer D'
%                type == 'dprimediff', D_{k}' - D_{ik}', for this routine k must be real
%                type == 'c', combined layer kernel coef(1) D + coef(2) S
%
%   varargin{1} - coef: length 2 array in the combined layer 
%                 formula, 2x2 matrix for all kernels
%                 otherwise does nothing
%
% Output:
%   submat - the evaluation of the selected kernel for the
%            provided sources and targets. the number of
%            rows equals the number of targets and the
%            number of columns equals the number of sources  
%
%

src = srcinfo.r(:,:); 
targ = targinfo.r(:,:);

[~, ns] = size(src);
[~, nt] = size(targ);

if (numel(zks) == 1)
    zk = zks(1);
else
    zk1 = zks(1);
    zk2 = zks(2);
end

switch lower(type)

case {'d', 'double'}
    srcnorm = srcinfo.n(:,:);
    [~, grad] = chnk.axissymhelm2d.green(zk, src, targ, origin);
    nx = repmat(srcnorm(1,:), nt, 1);
    ny = repmat(srcnorm(2,:), nt, 1);
    % Due to lack of translation invariance in r, no sign flip needed, 
    % as gradient is computed with repsect to r'
    submat = (grad(:,:,2).*nx - grad(:,:,3).*ny);
    fker = @(x, s, t, rns) fdlp(x, zk, s, t, rns, origin);
    for j=1:ns
        for i=1:nt
            rt = targ(1,i) + origin(1);
            dr = (src(1,j) - targ(1,i));
            dz = (src(2,j) - targ(2,i));
            r0   = sqrt(rt^2+(rt+dr)^2+dz^2);
            alph = (dr^2+dz^2)/r0^2;
            if alph > 2e-4 && alph < 0.2
                [x0, w0] = get_grid(zk, rt, dr, dz);
                fvals = fker(x0, src(:, j), targ(:,i), srcnorm(:,j));
                submat(i,j) = 2*w0.'*fvals;
            end
        end
    end 

case {'sp', 'sprime'}
    targnorm = targinfo.n(:,:);
    [~, grad] = chnk.axissymhelm2d.green(zk, src, targ, origin);

    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);
    submat = (grad(:,:,1).*nx + grad(:,:,3).*ny);


    fker = @(x, s, t, rnt) fsprime(x, zk, s, t, rnt, origin);
    for j=1:ns
        for i=1:nt
            rt = targ(1,i) + origin(1);
            dr = (src(1,j) - targ(1,i));
            dz = (src(2,j) - targ(2,i));
            r0   = sqrt(rt^2+(rt+dr)^2+dz^2);
            alph = (dr^2+dz^2)/r0^2;
            if alph > 2e-4 && alph < 0.2
                [x0, w0] = get_grid(zk, rt, dr, dz);
                fvals = fker(x0, src(:, j), targ(:,i), targnorm(:,i));
                submat(i,j) = 2*w0.'*fvals;
            end
        end
    end

case {'s', 'single'}
    submat = chnk.axissymhelm2d.green(zk, src, targ, origin);
    fker = @(x, s, t) fslp(x, zk, s, t, origin);
    for j=1:ns
        for i=1:nt
            rt = targ(1,i) + origin(1);
            dr = (src(1,j) - targ(1,i));
            dz = (src(2,j) - targ(2,i));
            r0   = sqrt(rt^2+(rt+dr)^2+dz^2);
            
            alph = (dr^2 + dz^2)/r0^2;
            
            if alph > 2e-4 && alph < 0.2
                [x0, w0] = get_grid(zk, rt, dr, dz);
                fvals = fker(x0, src(:, j), targ(:,i));
                submat(i,j) = 2*w0.'*fvals;
            end
        end
    end
    
case {'sdiff', 's_diff'}
    ifdiff = 1;
    submat = chnk.axissymhelm2d.green(zk, src, targ, origin, ifdiff);

case {'c', 'combined'}
    srcnorm = srcinfo.n(:,:); 
    coef = ones(2,1);
    if (nargin == 6); coef = varargin{1}; end
    nx = repmat(srcnorm(1,:), nt, 1);
    ny = repmat(srcnorm(2,:), nt, 1);
    [submats, grad] = chnk.axissymhelm2d.green(zk, src, targ, origin);
    % Due to lack of translation invariance in r, no sign flip needed, 
    % as gradient is computed with respect to r'
    submat = coef(1)*(grad(:,:,2).*nx - grad(:,:,3).*ny) + coef(2)*submats;
   
    fker = @(x, s, t, rns) coef(1)*fdlp(x, zk, s, t, rns, origin) + ...
           coef(2)*fslp(x, zk, s, t, origin);
    for j=1:ns
        for i=1:nt
            rt = targ(1,i) + origin(1);
            dr = (src(1,j) - targ(1,i));
            dz = (src(2,j) - targ(2,i));
            r0   = sqrt(rt^2+(rt+dr)^2+dz^2);
            alph = (dr^2+dz^2)/r0^2;
            if alph > 2e-4 && alph < 0.2
                [x0, w0] = get_grid(zk, rt, dr, dz);
                fvals = fker(x0, src(:, j), targ(:,i), srcnorm(:,j));
                submat(i,j) = 2*w0.'*fvals;
            end
        end
    end

case {'dp', 'dprime'}
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);
  [~,~,hess] = chnk.axissymhelm2d.green(zk, src, targ, origin);
  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);
  submat = hess(:,:,4).*nxsrc.*nxtarg - hess(:,:,5).*nysrc.*nxtarg ... 
              - hess(:,:,6).*nxsrc.*nytarg + hess(:,:,3).*nysrc.*nytarg;

case {'dp_re_diff', 'dprime_re_diff'}
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);
  ifdiff = 2;
  [~,~,hess] = chnk.axissymhelm2d.green(zk1, src, targ, origin,ifdiff);
  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);
  submat1 = hess(:,:,4).*nxsrc.*nxtarg - hess(:,:,5).*nysrc.*nxtarg ... 
              - hess(:,:,6).*nxsrc.*nytarg + hess(:,:,3).*nysrc.*nytarg;
  ifdiff = 2;
  [~,~,hess] = chnk.axissymhelm2d.green(zk2, src, targ, origin,ifdiff);
  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);
  submat2 = hess(:,:,4).*nxsrc.*nxtarg - hess(:,:,5).*nysrc.*nxtarg ... 
              - hess(:,:,6).*nxsrc.*nytarg + hess(:,:,3).*nysrc.*nytarg;

  submat = submat1-submat2;

case {'dpdiff', 'dp_diff', 'dprimediff', 'dprime_diff'}
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);
  ifdiff = 1;
  [~,~,hess] = chnk.axissymhelm2d.green(zk, src, targ, origin, ifdiff);
  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);
  submat = hess(:,:,4).*nxsrc.*nxtarg - hess(:,:,5).*nysrc.*nxtarg ...
      - hess(:,:,6).*nxsrc.*nytarg + hess(:,:,3).*nysrc.*nytarg;

    fker = @(x, s, t, rns, rnt) fdprimediff(x, zk, s, t, rns, rnt, origin);
    for j=1:ns
        for i=1:nt
            rt = targ(1,i) + origin(1);
            dr = (src(1,j) - targ(1,i));
            dz = (src(2,j) - targ(2,i));
            r0   = sqrt(rt^2+(rt+dr)^2+dz^2);
            alph = (dr^2+dz^2)/r0^2;
            if alph > 2e-4 && alph < 0.2
                [x0, w0] = get_grid(zk, rt, dr, dz);
                fvals = fker(x0, src(:, j), targ(:,i), srcnorm(:,j), ...
                    targnorm(:,i));
                submat(i,j) = 2*w0.'*fvals;
            end
        end
    end

case {'neu_rpcomb'}
  targnorm = targinfo.n(:,:);
  srcnorm = srcinfo.n(:,:);
  [~,gk,~,sikmat,gik,~,~,~,hessdiff] = ...
      chnk.axissymhelm2d.green_neu_all(zk, src, targ, origin);
  alpha = 1;
  if (nargin == 6); alpha = varargin{1}; end
  if (size(alpha) > 1)
    warning('Incorrect dimensions for coefs, using first component');
    alpha = alpha(1);
  end
    
  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);

  spmat = (gk(:,:,1).*nxtarg + gk(:,:,3).*nytarg);
  spikmat = (gik(:,:,1).*nxtarg + gik(:,:,3).*nytarg);
  dkdiffmat = hessdiff(:,:,4).*nxsrc.*nxtarg ...
      - hessdiff(:,:,5).*nysrc.*nxtarg ...
      - hessdiff(:,:,6).*nxsrc.*nytarg + hessdiff(:,:,3).*nysrc.*nytarg;

  rt = repmat(targ(1,:).',1,ns); rt = rt(:);
  rs = repmat(src(1,:),nt,1); rs = rs(:);
  dr = rs - rt;
  rt = rt + origin(1);
  dz = repmat(src(2,:),nt,1)-repmat(targ(2,:).',1,ns); dz = dz(:);
  r0 = rt.^2 + (rt + dr).^2 + dz.^2;
  alphs = (dr.^2 + dz.^2)./r0;
  iflag = (alphs > 2e-4).*(alphs < 0.2);
  if sum(iflag)

      % Assumes r,z are specified in meters, and, k is appropriately scaled
      rtmax = max(targ(1,:));
      rsmax = max(src(1,:));
      rmax = max(rtmax,rsmax);
      dr0 = 2e-4;
      dz0 = 2e-4;
      ppw = 10;
      [x0, w0] = get_grid(zk, rmax, dr0, dz0, ppw);
      iflag = reshape(iflag, [nt, ns]);
      sxhalf = sin(x0/2);
      sxhalf2 = sxhalf.*sxhalf;
      cx = 1- 2*sxhalf2;
      for j=1:ns
            for i=1:nt
                if iflag(i,j)
                    [fkp, fik, fikp, fkdiff] = get_neu_kers(zk, cx, sxhalf2, ...
                              src(:,j), targ(:,i), srcnorm(:,j), targnorm(:,i), origin);
                    spmat(i,j) = 2.*w0'*fkp;
                    sikmat(i,j) = 2.*w0'*fik;
                    spikmat(i,j) = 2.*w0'*fikp;
                    dkdiffmat(i,j) = 2.*w0'*fkdiff; 
                end
            end
      end
  end
  submat = zeros(3*nt, 3*ns);
  c1 = -1.0/(0.5 + 0.25*1i*alpha);
  c2 = 1i*alpha*c1;
  submat(1:3:end, 1:3:end) = c1*spmat;
  submat(1:3:end, 2:3:end) = c2*dkdiffmat;
  submat(1:3:end, 3:3:end) = c2*spikmat;
  submat(2:3:end, 1:3:end) = -sikmat;
  submat(3:3:end, 1:3:end) = -spikmat;
otherwise
    error('Unknown axissymmetric Helmholtz kernel type ''%s''.', type); 
end

end



function f = fslp (x, zk, s, t, o)
    rs = s(1); zs = s(2);
    rt = t(1); zt = t(2);
    
    r = sqrt((rs-rt).^2 + (zs-zt).^2 + 4*(rs+o(1)).*(rt+o(1)).*sin(x/2).^2);
    f = exp(1j*zk*r)/4/pi./r.*(rs + o(1));
end



function f = fdlp (x, zk, s, t, rns, o)
    rs = s(1); zs = s(2);
    rt = t(1); zt = t(2);
    
    rnd = ((rt +o(1)).*cos(x)  - (rs + o(1))).*rns(1) + (zt - zs).*rns(2);
    
    r = sqrt((rs-rt).^2 + (zs-zt).^2 + 4*(rs+o(1)).*(rt+o(1)).*sin(x/2).^2);
    f = rnd.*(1 - 1j*zk*r).*exp(1j*zk*r)/4/pi./r.^3.*(rs + o(1));
end



function f = fsprime (x, zk, s, t, rnt, o)
    rs = s(1); zs = s(2);
    rt = t(1); zt = t(2);
    
    rnd = ((rt + o(1)) - (rs + o(1)).*cos(x)).*rnt(1) + (zt - zs).*rnt(2);   
    r = sqrt((rs-rt).^2 + (zs-zt).^2 + 4*(rs + o(1)).*(rt + o(1)).*sin(x/2).^2);
    f = -rnd.*(1-1j*zk*r).*exp(1j*zk*r)/4/pi./r.^3.*(rs + o(1));
end



function f = fdprime (x, zk, s, t, rns, rnt, o)
    rs = s(1); zs = s(2);
    rt = t(1); zt = t(2);

    sxhalf = sin(x/2);
    sxhalf2 = sxhalf.*sxhalf;
    cx = 1-2*sxhalf2;
    
    rndt = ((rt + o(1)) - (rs + o(1)).*cx).*rnt(1) + (zt - zs).*rnt(2);
    rnds = ((rt + o(1)).*cx  - (rs + o(1))).*rns(1) + (zt - zs).*rns(2);
    rnsnt = rns(1)*rnt(1).*cx + rns(2)*rnt(2);
    
    r = sqrt((rs-rt).^2 + (zs-zt).^2 + 4*(rs+o(1)).*(rt+o(1)).*sxhalf2);
    f = -(rnsnt.*(1j*zk.*r-1)./r.^3 + ...
           rndt.*rnds.*(-zk^2.*r.^2 - 3*1j*zk.*r + 3)./r.^5).*exp(1j*zk*r)/4/pi.*(rs + o(1));
end

function [fkp, fik, fikp, fkdiff] = get_neu_kers(zk, cx, sxhalf2, s, t, rns, rnt, o)

    rs = s(1); zs = s(2);
    rt = t(1); zt = t(2);
    
        
    rndt = ((rt + o(1)) - (rs + o(1)).*cx).*rnt(1) + (zt - zs).*rnt(2);
    rnds = ((rt + o(1)).*cx  - (rs + o(1))).*rns(1) + (zt - zs).*rns(2);
    rnsnt = rns(1)*rnt(1).*cx + rns(2)*rnt(2);
    
    r = sqrt((rs-rt).^2 + (zs-zt).^2 + 4*(rs+o(1)).*(rt+o(1)).*sxhalf2);
    rinv = 1.0./r;
    rinv2 = rinv.*rinv;
    rinv3 = rinv.*rinv2;
    rinv5 = rinv3.*rinv2;
    
    afac = 1/4/pi.*(rs + o(1));
    efac = exp(1j*zk*r).*afac;
    efac_i = exp(-zk*r).*afac;
    
    fkp = -rndt.*(1-1j*zk*r).*rinv3.*efac;
    
    fik = efac_i.*rinv;
    fikp = -rndt.*(1 + zk*r).*rinv3.*efac_i;
    
    
    fkdiff = -(rnsnt.*(1j*zk.*r-1).*rinv3 + ...
           rndt.*rnds.*(-zk^2.*r.^2 - 3*1j*zk.*r + 3).*rinv5).*efac;
    fkdiff = fkdiff + (rnsnt.*(-zk.*r-1).*rinv3 + ...
           rndt.*rnds.*(zk^2.*r.^2 + 3*zk.*r + 3).*rinv5).*efac_i;
end


function f = fsdiff (x, zk, s, t, o)
    f = fslp(x, zk, s, t, o) - fslp(x, 1j*zk, s, t, o);
end



function f = fdprimediff (x, zk, s, t, rns, rnt, o)
    f1 = fdprime(x, zk, s, t, rns, rnt, o); 
    f2 = fdprime(x, 1j*zk, s, t, rns, rnt, o);
    f = f1 - f2;
end


function [xlegs, wlegs] = get_grid(zk, rt, dr, dz, ppw)

    if (nargin <= 4); ppw = 10; end
    persistent xlegloc wlegloc
    k = 16;
    if isempty(xlegloc) && isempty(wlegloc)
        [xlegloc,wlegloc] = lege.exps(k);
    end

    rs = rt + dr;
    drdiff = sqrt((rt+rs)^2  + dz^2) - sqrt(dr^2 + dz^2);
    npan = max(ceil(abs(zk)*drdiff*ppw/2/pi/k),3);
    
    h = pi/npan;
    tends = h:h:pi;
    
    dr0 = sqrt(dr^2 + dz^2);
    nref = max(ceil(-log(dr0/h)/log(2)),2);
    tends = [2.^(-nref:-1)*h tends];
    tstarts = [0 tends(1:end-1)];
    
    xlegs = tstarts + (xlegloc+1)/2*(tends-tstarts);
    wlegs = wlegloc/2*(tends-tstarts);
    xlegs = xlegs(:);
    wlegs = wlegs(:);
end

