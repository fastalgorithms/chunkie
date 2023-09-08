function submat = kern(zk, srcinfo, targinfo, origin, type, varargin)
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
%                type == 'ddiff', D_{k}' - D_{ik}', for this routine k must be real
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

src = srcinfo.r; 
targ = targinfo.r;

[~, ns] = size(src);
[~, nt] = size(targ);

if strcmpi(type, 'd')
    srcnorm = srcinfo.n;
    [~, grad] = chnk.axissymhelm2d.green(zk, src, targ, origin);
    nx = repmat(srcnorm(1,:), nt, 1);
    ny = repmat(srcnorm(2,:), nt, 1);
    % Due to lack of translation invariance in r, no sign flip needed, 
    % as gradient is computed with repsect to r'
    submat = (grad(:,:,2).*nx - grad(:,:,3).*ny);
    fker = @(x, s, t, rns) fdlp(x, zk, s, t, rns);
    for j=1:ns
        for i=1:nt
            rt = targ(1,i);
            dr = (src(1,j) - targ(1,i));
            dz = (src(2,j) - targ(2,i));
            r0   = sqrt(rt^2+(rt+dr)^2+dz^2);
            alph = (dr^2+dz^2)/r0^2;
            if alph > 2e-4
                [x0, w0] = get_grid(zk, rt, dr, dz);
                fvals = fker(x0, src(:, j), targ(:,i), srcnorm(:,j));
                submat(i,j) = 2*w0.'*fvals;
                % submat(i,j) = integral(@(x) fker(x, src(:,j), ...
                %       targ(:,i), srcnorm(:,j)), 0, 2*pi, ...
                %       'AbsTol',1e-14,'RelTol',1e-10); 
            end
        end
    end 
end

if strcmpi(type, 'sprime')
    targnorm = targinfo.n;
    [~, grad] = chnk.axissymhelm2d.green(zk, src, targ, origin);

    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);
    submat = (grad(:,:,1).*nx + grad(:,:,3).*ny);


    fker = @(x, s, t, rnt) fsprime(x, zk, s, t, rnt);
    for j=1:ns
        for i=1:nt
            rt = targ(1,i);
            dr = (src(1,j) - targ(1,i));
            dz = (src(2,j) - targ(2,i));
            r0   = sqrt(rt^2+(rt+dr)^2+dz^2);
            alph = (dr^2+dz^2)/r0^2;
            if alph > 2e-4
                [x0, w0] = get_grid(zk, rt, dr, dz);
                fvals = fker(x0, src(:, j), targ(:,i), targnorm(:,i));
                submat(i,j) = 2*w0.'*fvals;
                % submat(i,j) = integral(@(x) fker(x, src(:,j), ...
                %       targ(:,i), targnorm(:,i)), 0, 2*pi, ...
                %       'AbsTol',1e-14,'RelTol',1e-10); 
            end
        end
    end

end

if strcmpi(type, 's')
    submat = chnk.axissymhelm2d.green(zk, src, targ, origin);
    fker = @(x, s, t) fslp(x, zk, s, t);
    for j=1:ns
        for i=1:nt
            rt = targ(1,i);
            dr = (src(1,j) - targ(1,i));
            dz = (src(2,j) - targ(2,i));
            r0   = sqrt(rt^2+(rt+dr)^2+dz^2);
            
            alph = (dr^2 + dz^2)/r0^2;
            
            if alph > 2e-4
                [x0, w0] = get_grid(zk, rt, dr, dz);
                fvals = fker(x0, src(:, j), targ(:,i));
                submat(i,j) = 2*w0.'*fvals;
                % submat(i,j) = integral(@(x) fker(x, src(:,j), ...
                %       targ(:,i)), 0, 2*pi, ...
                %       'AbsTol',1e-14,'RelTol',1e-10); 
            end
        end
    end
    
end


if strcmpi(type, 'sdiff')
    ifdiff = 1;
    submat = chnk.axissymhelm2d.green(zk, src, targ, origin, ifdiff);
end

if strcmpi(type, 'c')
    srcnorm = srcinfo.n; 
    coef = ones(2,1);
    if (nargin == 6); coef = varargin{1}; end
    nx = repmat(srcnorm(1,:), nt, 1);
    ny = repmat(srcnorm(2,:), nt, 1);
    [submats, grad] = chnk.axissymhelm2d.green(zk, src, targ, origin);
    % Due to lack of translation invariance in r, no sign flip needed, 
    % as gradient is computed with respect to r'
    submat = coef(1)*(grad(:,:,2).*nx - grad(:,:,3).*ny) + coef(2)*submats;
   
    fker = @(x, s, t, rns) coef(1)*fdlp(x, zk, s, t, rns) + coef(2)*fslp(x, zk, s, t);
    for j=1:ns
        for i=1:nt
            rt = targ(1,i);
            dr = (src(1,j) - targ(1,i));
            dz = (src(2,j) - targ(2,i));
            r0   = sqrt(rt^2+(rt+dr)^2+dz^2);
            alph = (dr^2+dz^2)/r0^2;
            if alph > 2e-4
                [x0, w0] = get_grid(zk, rt, dr, dz);
                fvals = fker(x0, src(:, j), targ(:,i), srcnorm(:,j));
                submat(i,j) = 2*w0.'*fvals;
                % submat(i,j) = integral(@(x) fker(x, src(:,j), ...
                %       targ(:,i), srcnorm(:,j)), 0, 2*pi, ...
                %       'AbsTol',1e-14,'RelTol',1e-10); 
            end
        end
    end
end


if strcmpi(type, 'dprime')
  targnorm = targinfo.n;
  srcnorm = srcinfo.n;
  [~,~,hess] = chnk.axissymhelm2d.green(zk, src, targ, origin);
  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);
  submat = hess(:,:,4).*nxsrc.*nxtarg - hess(:,:,5).*nysrc.*nxtarg ... 
              - hess(:,:,6).*nxsrc.*nytarg + hess(:,:,3).*nysrc.*nytarg;
end



if strcmpi(type, 'dprimediff')
  targnorm = targinfo.n;
  srcnorm = srcinfo.n;
  ifdiff = 1;
  [~,~,hess] = chnk.axissymhelm2d.green(zk, src, targ, origin, ifdiff);
  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);
  submat = hess(:,:,4).*nxsrc.*nxtarg - hess(:,:,5).*nysrc.*nxtarg ...
      - hess(:,:,6).*nxsrc.*nytarg + hess(:,:,3).*nysrc.*nytarg;

    fker = @(x, s, t, rns, rnt) fdprimediff(x, zk, s, t, rns, rnt);
    for j=1:ns
        for i=1:nt
            rt = targ(1,i);
            dr = (src(1,j) - targ(1,i));
            dz = (src(2,j) - targ(2,i));
            r0   = sqrt(rt^2+(rt+dr)^2+dz^2);
            alph = (dr^2+dz^2)/r0^2;
            if alph > 2e-4
                [x0, w0] = get_grid(zk, rt, dr, dz);
                fvals = fker(x0, src(:, j), targ(:,i), srcnorm(:,j), ...
                    targnorm(:,i));
                submat(i,j) = 2*w0.'*fvals;
                % submat(i,j) = integral(@(x) fker(x, src(:,j), ...
                %       targ(:,i), srcnorm(:,j), targnorm(:,i)), 0, 2*pi, ...
                %       'AbsTol',1e-14,'RelTol',1e-10); 
            end
        end
    end

end

end



function f = fslp (x, zk, s, t)
    rs = s(1); zs = s(2);
    rt = t(1); zt = t(2);
    
    % r = sqrt(rs.^2 + rt.^2 - 2*rs.*rt.*cos(x) + (zs-zt).^2);
    r = sqrt((rs-rt).^2 + (zs-zt).^2 + 4*rs.*rt.*sin(x/2).^2);
    f = exp(1j*zk*r)/4/pi./r.*rs;
end



function f = fdlp (x, zk, s, t, rns)
    rs = s(1); zs = s(2);
    rt = t(1); zt = t(2);
    
    rnd = (rt.*cos(x)  - rs).*rns(1) + (zt - zs).*rns(2);
    
    % r = sqrt(rs.^2 + rt.^2 - 2*rs.*rt.*cos(x) + (zs-zt).^2);
    r = sqrt((rs-rt).^2 + (zs-zt).^2 + 4*rs.*rt.*sin(x/2).^2);
    f = rnd.*(1 - 1j*zk*r).*exp(1j*zk*r)/4/pi./r.^3.*rs;
end



function f = fsprime (x, zk, s, t, rnt)
    rs = s(1); zs = s(2);
    rt = t(1); zt = t(2);
    
    rnd = (rt - rs.*cos(x)).*rnt(1) + (zt - zs).*rnt(2);
    
    % r = sqrt(rs.^2 + rt.^2 - 2*rs.*rt.*cos(x) + (zs-zt).^2);
    r = sqrt((rs-rt).^2 + (zs-zt).^2 + 4*rs.*rt.*sin(x/2).^2);
    f = -rnd.*(1-1j*zk*r).*exp(1j*zk*r)/4/pi./r.^3.*rs;
end



function f = fdprime (x, zk, s, t, rns, rnt)
    rs = s(1); zs = s(2);
    rt = t(1); zt = t(2);
    
    rndt = (rt - rs.*cos(x)).*rnt(1) + (zt - zs).*rnt(2);
    rnds = (rt.*cos(x) - rs).*rns(1) + (zt - zs).*rns(2);
    rnsnt = rns(1)*rnt(1).*cos(x) + rns(2)*rnt(2);
    
    % r = sqrt(rs.^2 + rt.^2 - 2*rs.*rt.*cos(x) + (zs-zt).^2);
    r = sqrt((rs-rt).^2 + (zs-zt).^2 + 4*rs.*rt.*sin(x/2).^2);
    f = -(rnsnt.*(1j*zk.*r-1).*exp(1j*zk*r)/4/pi./r.^3 + ...
           rndt.*rnds.*(-zk^2.*r.^2 - 3*1j*zk.*r + 3).*exp(1j*zk*r)/4/pi./r.^5).*rs;
end


function f = fsdiff (x, zk, s, t)
    f = fslp(x, zk, s, t) - fslp(x, 1j*zk, s, t);
end



function f = fdprimediff (x, zk, s, t, rns, rnt)
    f1 = fdprime(x, zk, s, t, rns, rnt); 
    f2 = fdprime(x, 1j*zk, s, t, rns, rnt);
    f = f1 - f2;
end


function [xlegs, wlegs] = get_grid(zk, rt, dr, dz, ppw)

    if (nargin <= 4); ppw = 30; end
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
    nref = nref + 2;
    if (nref>20); fprintf('nref = %d   dr0 = %d\n',nref,dr0); end
    tends = [2.^(-nref:-1)*h tends];
    tstarts = [0 tends(1:end-1)];
    
    xlegs = tstarts + (xlegloc+1)/2*(tends-tstarts);
    wlegs = wlegloc/2*(tends-tstarts);
    xlegs = xlegs(:);
    wlegs = wlegs(:);
end

