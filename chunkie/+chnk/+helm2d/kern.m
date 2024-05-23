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
% D'(x,y) = \nabla_{n_x} \nabla_{n_y} G(x,y)
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
%                type == 'dprime', normal derivative of double layer D'
%                type == 'c', combined layer kernel coef(1) D + coef(2) S
%                type == 'stau', tangential derivative of single layer
%                type == 'all', returns all four layer potentials, 
%                       [coef(1,1)*D coef(1,2)*S; coef(2,1)*D' coef(2,2)*S']
%                type == 'c2trans' returns the combined field, and the 
%                          normal derivative of the combined field
%                        [coef(1)*D + coef(2)*S; coef(1)*D' + coef(2)*S']
%                type == 'trans_rep' returns the potential corresponding
%                           to the transmission representation
%                        [coef(1)*D coef(2)*S]
%                type == 'trans_rep_prime' returns the normal derivative
%                          corresponding to the transmission representation
%                        [coef(1)*D' coef(2)*S']
%                type == 'trans_rep_grad' returns the gradient corresponding
%                         to the transmission representation
%                        [coef(1)*d_x D coef(2)*d_x S;
%                         coef(1)*d_y D coef(2)*d_y S]
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
% see also CHNK.HELM2D.GREEN
  
src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);

% double layer
if strcmpi(type,'d')
  srcnorm = srcinfo.n;
  [~,grad] = chnk.helm2d.green(zk,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  submat = -(grad(:,:,1).*nx + grad(:,:,2).*ny);
end

% normal derivative of single layer
if strcmpi(type,'sprime')
  targnorm = targinfo.n;
  [~,grad] = chnk.helm2d.green(zk,src,targ);
  nx = repmat((targnorm(1,:)).',1,ns);
  ny = repmat((targnorm(2,:)).',1,ns);
  submat = (grad(:,:,1).*nx + grad(:,:,2).*ny);
end


% Tangential derivative of single layer
if strcmpi(type,'stau')
  targtan = targinfo.d;
  [~,grad] = chnk.helm2d.green(zk,src,targ);
  dx = repmat((targtan(1,:)).',1,ns);
  dy = repmat((targtan(2,:)).',1,ns);
  ds = sqrt(dx.*dx+dy.*dy);
  submat = (grad(:,:,1).*dx + grad(:,:,2).*dy)./ds;
end

% single layer
if strcmpi(type,'s')
  submat = chnk.helm2d.green(zk,src,targ);
end

% normal derivative of double layer
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

% Combined field 
if strcmpi(type,'c')
  srcnorm = srcinfo.n;
  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end
  [submats,grad] = chnk.helm2d.green(zk,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  submatd = -(grad(:,:,1).*nx + grad(:,:,2).*ny);
  submat = coef(1)*submatd + coef(2)*submats;
end

% normal derivative of combined field
if strcmpi(type,'cprime')
  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end;
  targnorm = targinfo.n;
  srcnorm = srcinfo.n;
  
  submat = zeros(nt,ns);

  % Get gradient and hessian info
  [submats,grad,hess] = chnk.helm2d.green(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);

  % D'
  submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
      + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
  % S'
  submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);

  submat = coef(1)*submatdp + coef(2)*submatsp;
end


% Dirichlet and neumann data corresponding to combined field
if strcmpi(type,'c2trans') 
  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end;
  targnorm = targinfo.n;
  srcnorm = srcinfo.n;

  % Get gradient and hessian info
  [submats,grad,hess] = chnk.helm2d.green(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);

  % D'
  submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
      + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
  % S'
  submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
  % D
  submatd = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);

  submat = zeros(2*nt,ns);

  submat(1:2:2*nt,:) = coef(1)*submatd + coef(2)*submats;
  submat(2:2:2*nt,:) = coef(1)*submatdp + coef(2)*submatsp;

end


% all kernels, [c11 D, c12 S; c21 D', c22 S'] 
if strcmpi(type,'all')
 
  targnorm = targinfo.n;
  srcnorm = srcinfo.n;
  cc = varargin{1};
  
  submat = zeros(2*nt,2*ns);

  % Get gradient and hessian info
  [submats,grad,hess] = chnk.helm2d.green(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);

  submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
      + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
  % S'
  submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
  % D
  submatd  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);
  
  
  submat(1:2:2*nt,1:2:2*ns) = submatd*cc(1,1);
  submat(1:2:2*nt,2:2:2*ns) = submats*cc(1,2);
  submat(2:2:2*nt,1:2:2*ns) = submatdp*cc(2,1);
  submat(2:2:2*nt,2:2:2*ns) = submatsp*cc(2,2);
end


% clamped plate kernel for modified biharmonic problem
if strcmpi(type, 'clamped-plate')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;

   [hess, third, forth] = chnk.helm2d.thirdforth_derivatives(zk, src, targ);           % Hankel part
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   
   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   zkimag = (1i)*zk;
   [hessK, thirdK, forthK] = chnk.helm2d.thirdforth_derivatives(zkimag, src, targ);     % modified bessel K part

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);
    
    ds = sqrt(dx.*dx+dy.*dy);

    taux = dx./ds;
    tauy = dy./ds;

    Kxx = 1/(2*zk^2).*(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
       third(:, :, 4).*(3*nx.*ny.*ny) + third(:, :, 6).*(ny.*ny.*ny)) - ...
        1/(2*zk^2).*(thirdK(:, :, 1).*(nx.*nx.*nx) + thirdK(:, :, 2).*(3*nx.*nx.*ny) +...
       thirdK(:, :, 4).*(3*nx.*ny.*ny) + thirdK(:, :, 6).*(ny.*ny.*ny)) + ...
       (3/(2*zk^2).*(third(:, :, 1).*(nx.*taux.*taux) + third(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       third(:, :, 4).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + third(:, :, 6).*(ny.*tauy.*tauy))-...
       3/(2*zk^2).*(thirdK(:, :, 1).*(nx.*taux.*taux) + thirdK(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       thirdK(:, :, 4).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + thirdK(:, :, 6).*(ny.*tauy.*tauy)));  % G_{ny ny ny} + 3G_{ny tauy tauy}

   Kxy = -(1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny))-...
           1/(2*zk^2).*(hessK(:, :, 1).*(nx.*nx) + hessK(:, :, 2).*(2*nx.*ny) + hessK(:, :, 3).*(ny.*ny)))+...
          (1/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy))-...
           1/(2*zk^2).*(hessK(:, :, 1).*(taux.*taux) + hessK(:, :, 2).*(2*taux.*tauy) + hessK(:, :, 3).*(tauy.*tauy))); % -G_{ny ny}  + G_{tauy tauy}

   Kyx = 1/(2*zk^2).*(forth(:, :, 1).*(nx.*nx.*nx.*nxtarg) + forth(:, :, 2).*(nx.*nx.*nx.*nytarg + 3*nx.*nx.*ny.*nxtarg) + ...
          forth(:, :, 4).*(3*nx.*nx.*ny.*nytarg + 3*nx.*ny.*ny.*nxtarg) + forth(:, :, 6).*(3*nx.*ny.*ny.*nytarg +ny.*ny.*ny.*nxtarg)+...
          forth(:, :, 8).*(ny.*ny.*ny.*nytarg)) - 1/(2*zk^2).*(forthK(:, :, 1).*(nx.*nx.*nx.*nxtarg) + forthK(:, :, 2).*(nx.*nx.*nx.*nytarg + 3*nx.*nx.*ny.*nxtarg) + ...
          forthK(:, :, 4).*(3*nx.*nx.*ny.*nytarg + 3*nx.*ny.*ny.*nxtarg) + forthK(:, :, 6).*(3*nx.*ny.*ny.*nytarg +ny.*ny.*ny.*nxtarg)+...
          forthK(:, :, 8).*(ny.*ny.*ny.*nytarg)) + ...
          (3/(2*zk^2).*(forth(:, :, 1).*(nx.*taux.*taux.*nxtarg)+ forth(:, :, 2).*(nx.*taux.*taux.*nytarg + 2*nx.*taux.*tauy.*nxtarg + ny.*taux.*taux.*nxtarg) +...
          forth(:, :, 4).*(2*nx.*taux.*tauy.*nytarg + ny.*taux.*taux.*nytarg + nx.*tauy.*tauy.*nxtarg + 2*ny.*taux.*tauy.*nxtarg) + ...
          forth(:, :, 6).*(nx.*tauy.*tauy.*nytarg +2*ny.*taux.*tauy.*nytarg + ny.*tauy.*tauy.*nxtarg) + forth(:, :, 8).*(ny.*tauy.*tauy.*nytarg))-...
          3/(2*zk^2).*(forthK(:, :, 1).*(nx.*taux.*taux.*nxtarg)+ forthK(:, :, 2).*(nx.*taux.*taux.*nytarg + 2*nx.*taux.*tauy.*nxtarg + ny.*taux.*taux.*nxtarg) +...
          forthK(:, :, 4).*(2*nx.*taux.*tauy.*nytarg + ny.*taux.*taux.*nytarg + nx.*tauy.*tauy.*nxtarg + 2*ny.*taux.*tauy.*nxtarg) + ...
          forthK(:, :, 6).*(nx.*tauy.*tauy.*nytarg +2*ny.*taux.*tauy.*nytarg + ny.*tauy.*tauy.*nxtarg) + forthK(:, :, 8).*(ny.*tauy.*tauy.*nytarg)));

   Kyy = 1/(2*zk^2).*(third(:,:, 1).*(nx.*nx.*nxtarg) +third(:, :, 2).*(nx.*nx.*nytarg + 2*nx.*ny.*nxtarg) + third(:, :, 4).*(2*nx.*ny.*nytarg + ny.*ny.*nxtarg)+...
         third(:, :,6).*(ny.*ny.*nytarg)) - 1/(2*zk^2).*(thirdK(:,:, 1).*(nx.*nx.*nxtarg) +thirdK(:, :, 2).*(nx.*nx.*nytarg + 2*nx.*ny.*nxtarg) + thirdK(:, :, 4).*(2*nx.*ny.*nytarg + ny.*ny.*nxtarg)+...
         thirdK(:, :,6).*(ny.*ny.*nytarg)) - 1/(2*zk^2).*(third(:,:, 1).*(taux.*taux.*nxtarg) +third(:, :, 2).*(taux.*taux.*nytarg + 2*taux.*tauy.*nxtarg) + third(:, :, 4).*(2*taux.*tauy.*nytarg + tauy.*tauy.*nxtarg)+...
         third(:, :,6).*(tauy.*tauy.*nytarg)) + 1/(2*zk^2).*(thirdK(:,:, 1).*(taux.*taux.*nxtarg) +thirdK(:, :, 2).*(taux.*taux.*nytarg + 2*taux.*tauy.*nxtarg) + thirdK(:, :, 4).*(2*taux.*tauy.*nytarg + tauy.*tauy.*nxtarg)+...
         thirdK(:, :,6).*(tauy.*tauy.*nytarg));

  submat = zeros(2*nt,2*ns);
  
  submat(1:2:end,1:2:end) = Kxx;
  submat(1:2:end,2:2:end) = Kxy;
    
  submat(2:2:end,1:2:end) = Kyx;
  submat(2:2:end,2:2:end) = Kyy;
end

% free plate kernel for the modified biharmonic problem (K11 with no
% hilbert transform subtraction, K12 kernel, K22 with no curvature part) K21 is
% handled in a separate type.
if strcmpi(type, 'free plate first part')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};
   



   [hess, third, ~] = chnk.helm2d.thirdforth_derivatives(zk, src, targ);            % Hankel part
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   zkimag = (1i)*zk;
   [hessK, thirdK, ~] = chnk.helm2d.thirdforth_derivatives(zkimag, src, targ);     % modified bessel K part
   

  
   

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   
   
   
   K11 = 1/(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nx) + third(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
        third(:, :, 4).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
        third(:, :, 6).*(nytarg.*nytarg.*ny)) - ...
        1/(2*zk^2).*(thirdK(:, :, 1).*(nxtarg.*nxtarg.*nx) + thirdK(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
        thirdK(:, :, 4).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
        thirdK(:, :, 6).*(nytarg.*nytarg.*ny)) + ...
       coefs(1)/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + third(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
        third(:, :, 4).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
        third(:, :, 6).*(tauytarg.*tauytarg.*ny)) - ...
        coefs(1)/(2*zk^2).*(thirdK(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + thirdK(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
        thirdK(:, :, 4).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
        thirdK(:, :, 6).*(tauytarg.*tauytarg.*ny));  % first kernel with no hilbert transforms (G_{nx nx ny + nu G_{taux taux ny}).
       
    
       

   K12 =  1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg)-...
           1/(2*zk^2).*(hessK(:, :, 1).*nxtarg.*nxtarg + hessK(:, :, 2).*(2*nxtarg.*nytarg) + hessK(:, :, 3).*nytarg.*nytarg)+...
           coefs(1)/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg)-...
           coefs(1)/(2*zk^2).*(hessK(:, :, 1).*tauxtarg.*tauxtarg + hessK(:, :, 2).*(2*tauxtarg.*tauytarg) + ...
           hessK(:, :, 3).*tauytarg.*tauytarg) ;    % G_{nx nx} + nu G_{taux taux}
   
   K21 = 0;




   K22 = -1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nxtarg) + third(:, :, 2).*(3*nxtarg.*nxtarg.*nytarg) +...
       third(:, :, 4).*(3*nxtarg.*nytarg.*nytarg) + third(:, :, 6).*(nytarg.*nytarg.*nytarg)) + ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(nxtarg.*nxtarg.*nxtarg) + thirdK(:, :, 2).*(3*nxtarg.*nxtarg.*nytarg)+...
        thirdK(:, :, 4).*(3*nxtarg.*nytarg.*nytarg) + thirdK(:, :, 6).*(nytarg.*nytarg.*nytarg)) -...
        (2-coefs(1))/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg) + third(:, :, 2).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
        third(:, :, 4).*(2*tauxtarg.*tauytarg.*nytarg + tauytarg.*tauytarg.*nxtarg) +...
        + third(:, :, 6).*(tauytarg.*tauytarg.*nytarg)) + ...
        (2-coefs(1))/(2*zk^2).*(thirdK(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg) + thirdK(:, :, 2).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
        thirdK(:, :, 4).*(2*tauxtarg.*tauytarg.*nytarg + tauytarg.*tauytarg.*nxtarg) +...
        + thirdK(:, :, 6).*(tauytarg.*tauytarg.*nytarg)); % G_{nx nx nx} + (2-nu) G_{taux taux nx}
    
  submat = zeros(2*nt,2*ns);
  
  submat(1:2:end,1:2:end) = K11;
  submat(1:2:end,2:2:end) = K12;
    
  submat(2:2:end,1:2:end) = K21;
  submat(2:2:end,2:2:end) = K22;

end


% free plate kernel K21 for the modified biharmonic problem. This part
% handles the singularity subtraction and swap the evaluation to its
% asymptotic expansions if the targets and sources are close. 

if strcmpi(type, 'free plate K21 first part')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};
   R = 0.001; % radius of kernel replacement



   [~, ~, forth] = chnk.helm2d.thirdforth_derivatives(zk, src, targ);            % Hankel part
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   zkimag = (1i)*zk;
   [~, ~, forthK] = chnk.helm2d.thirdforth_derivatives(zkimag, src, targ);     % modified bessel K part
    

  
   

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;
   dists = sqrt(r2);
   inds = and((dists < R),(dists > 0));

   submat  = 1/(2*zk^2).*(forth(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*nx) + forth(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*ny + 3*nxtarg.*nxtarg.*nytarg.*nx) + ...
          forth(:, :, 4).*(3*nxtarg.*nxtarg.*nytarg.*ny + 3*nxtarg.*nytarg.*nytarg.*nx) + forth(:, :, 6).*(3*nxtarg.*nytarg.*nytarg.*ny +nytarg.*nytarg.*nytarg.*nx)+...
          forth(:, :, 8).*(nytarg.*nytarg.*nytarg.*ny)) -...
          1/(2*zk^2).*(forthK(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*nx) + forthK(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*ny + 3*nxtarg.*nxtarg.*nytarg.*nx) + ...
          forthK(:, :, 4).*(3*nxtarg.*nxtarg.*nytarg.*ny + 3*nxtarg.*nytarg.*nytarg.*nx) + forthK(:, :, 6).*(3*nxtarg.*nytarg.*nytarg.*ny +nytarg.*nytarg.*nytarg.*nx)+...
          forthK(:, :, 8).*(nytarg.*nytarg.*nytarg.*ny)) + (3/(4*pi)).*(nxtarg.*nx + nytarg.*ny)./r2 + ...
          ((2-coefs(1))/(2*zk^2).*(forth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*nx) + forth(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*ny + tauxtarg.*tauxtarg.*nytarg.*nx + 2*tauxtarg.*tauytarg.*nxtarg.*nx)+...
          forth(:, :, 4).*(tauxtarg.*tauxtarg.*nytarg.*ny + 2*tauxtarg.*tauytarg.*nxtarg.*ny + tauytarg.*tauytarg.*nxtarg.*nx + 2*tauxtarg.*tauytarg.*nytarg.*nx)+...
          forth(:, :, 6).*(tauytarg.*tauytarg.*nxtarg.*ny + 2*tauxtarg.*tauytarg.*nytarg.*ny + tauytarg.*tauytarg.*nytarg.*nx) +...
          forth(:, :, 8).*(tauytarg.*tauytarg.*nytarg.*ny)) -...
          (2-coefs(1))/(2*zk^2).*(forthK(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*nx) + forthK(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*ny + tauxtarg.*tauxtarg.*nytarg.*nx + 2*tauxtarg.*tauytarg.*nxtarg.*nx)+...
          forthK(:, :, 4).*(tauxtarg.*tauxtarg.*nytarg.*ny + 2*tauxtarg.*tauytarg.*nxtarg.*ny + tauytarg.*tauytarg.*nxtarg.*nx + 2*tauxtarg.*tauytarg.*nytarg.*nx)+...
          forthK(:, :, 6).*(tauytarg.*tauytarg.*nxtarg.*ny + 2*tauxtarg.*tauytarg.*nytarg.*ny + tauytarg.*tauytarg.*nytarg.*nx) +...
          forthK(:, :, 8).*(tauytarg.*tauytarg.*nytarg.*ny)) +...
          (2-coefs(1))/(4*pi).*(nxtarg.*nx + nytarg.*ny)./r2 - ...
          ((2-coefs(1))/(2*pi)).*(nxtarg.*nx + nytarg.*ny).*(rx.*tauxtarg + ry.*tauytarg).*(rx.*tauxtarg + ry.*tauytarg)./(r2.^2)) - ...
          ((1+coefs(1))/(4*pi)).*((taux.*tauxtarg + tauy.*tauytarg)./(r2) - 2*(rx.*tauxtarg + ry.*tauytarg).*(rx.*taux + ry.*tauy)./(r2.^2))-...
          (3/(4*pi)).*(nxtarg.*nx + nytarg.*ny)./r2  - (2-coefs(1))/(4*pi).*(nxtarg.*nx + nytarg.*ny)./r2 + ...
          (2-coefs(1))/(2*pi).*(nxtarg.*nx + nytarg.*ny).*(rx.*tauxtarg + ry.*tauytarg).*(rx.*tauxtarg + ry.*tauytarg)./(r2.^2);
   % third kernel with singularity subtraction and H[delta']
 
   if any(inds(:))                                                      % if statement for evalauation points that are closed to the boundary
       temp = 1i*imag(submat(inds));
    
       nx = nx(inds); ny = ny(inds); 
      
       nxtarg = nxtarg(inds);
       nytarg = nytarg(inds);
    
       dx = dx(inds);
       dy = dy(inds);
    
       dx1 = dx1(inds);
       dy1 = dy1(inds);
    
       ds = sqrt(dx.*dx+dy.*dy);
       ds1 = sqrt(dx1.*dx1+dy1.*dy1); 
    
       taux = dx./ds;                                                                       % normalization
       tauy = dy./ds;
    
       tauxtarg = dx1./ds1;
       tauytarg = dy1./ds1;
    
       rx = rx(inds);
       ry = ry(inds);
       r2 = r2(inds);
    
       tauxtargnsrc = tauxtarg.*nx + tauytarg.*ny;
       tauxtargntarg = tauxtarg.*nxtarg + tauytarg.*nytarg;
    
       normaldot = nx.*nxtarg + ny.*nytarg;
       tangentdot = taux.*tauxtarg + tauy.*tauytarg;
   
       rn = rx.*nx + ry.*ny;

       rntarg = rx.*nxtarg + ry.*nytarg;

       rtau = rx.*taux + ry.*tauy;

       rtautarg = rx.*tauxtarg + ry.*tauytarg;
      
       

       [smoothfirst, smoothsecond] = chnk.helm2d.smooth_part_derivatives(zk, rx, ry);

       eulerconstantpart = (smoothfirst(:, :, 1)).*(nxtarg.*nxtarg.*nxtarg.*nx) +...
           (smoothfirst(:, :, 2)).*(nxtarg.*nxtarg.*nxtarg.*ny + 3*nxtarg.*nxtarg.*nytarg.*nx) + ...
          (smoothfirst(:, :, 3)).*(3*nxtarg.*nxtarg.*nytarg.*ny + 3*nxtarg.*nytarg.*nytarg.*nx) +...
          (smoothfirst(:, :, 4)).*(3*nxtarg.*nytarg.*nytarg.*ny +nytarg.*nytarg.*nytarg.*nx)+...
          (smoothfirst(:, :, 5)).*(nytarg.*nytarg.*nytarg.*ny)+ ...
          (2-coefs(1)).*((smoothfirst(:, :, 1)).*(tauxtarg.*tauxtarg.*nxtarg.*nx) +...
          (smoothfirst(:, :, 2)).*(tauxtarg.*tauxtarg.*nxtarg.*ny + tauxtarg.*tauxtarg.*nytarg.*nx + 2*tauxtarg.*tauytarg.*nxtarg.*nx)+...
          (smoothfirst(:, :, 3)).*(tauxtarg.*tauxtarg.*nytarg.*ny + 2*tauxtarg.*tauytarg.*nxtarg.*ny + tauytarg.*tauytarg.*nxtarg.*nx + 2*tauxtarg.*tauytarg.*nytarg.*nx)+...
          (smoothfirst(:, :, 4)).*(tauytarg.*tauytarg.*nxtarg.*ny + 2*tauxtarg.*tauytarg.*nytarg.*ny + tauytarg.*tauytarg.*nytarg.*nx) +...
          (smoothfirst(:, :, 5)).*(tauytarg.*tauytarg.*nytarg.*ny));

       puresmoothpart = (smoothsecond(:, :, 1)).*(nxtarg.*nxtarg.*nxtarg.*nx) +...
           (smoothsecond(:, :, 2)).*(nxtarg.*nxtarg.*nxtarg.*ny + 3*nxtarg.*nxtarg.*nytarg.*nx) + ...
          (smoothsecond(:, :, 3)).*(3*nxtarg.*nxtarg.*nytarg.*ny + 3*nxtarg.*nytarg.*nytarg.*nx) +...
          (smoothsecond(:, :, 4)).*(3*nxtarg.*nytarg.*nytarg.*ny +nytarg.*nytarg.*nytarg.*nx)+...
          (smoothsecond(:, :, 5)).*(nytarg.*nytarg.*nytarg.*ny)+ ...
          (2-coefs(1)).*((smoothsecond(:, :, 1)).*(tauxtarg.*tauxtarg.*nxtarg.*nx) +...
          (smoothsecond(:, :, 2)).*(tauxtarg.*tauxtarg.*nxtarg.*ny + tauxtarg.*tauxtarg.*nytarg.*nx + 2*tauxtarg.*tauytarg.*nxtarg.*nx)+...
          (smoothsecond(:, :, 3)).*(tauxtarg.*tauxtarg.*nytarg.*ny + 2*tauxtarg.*tauytarg.*nxtarg.*ny + tauytarg.*tauytarg.*nxtarg.*nx + 2*tauxtarg.*tauytarg.*nytarg.*nx)+...
          (smoothsecond(:, :, 4)).*(tauytarg.*tauytarg.*nxtarg.*ny + 2*tauxtarg.*tauytarg.*nytarg.*ny + tauytarg.*tauytarg.*nytarg.*nx) +...
          (smoothsecond(:, :, 5)).*(tauytarg.*tauytarg.*nytarg.*ny));

       submat(inds) = 3/(2*pi).*(rn.*rntarg./(r2.^2)) + 3/(2*pi).*(rntarg.^2)./(r2.^2) -...
            (2/pi).*(rn.*(rntarg.^3))./(r2.^3) + ...
            (2-coefs(1)).*(-1/(2*pi).*(tauxtargnsrc).*(tauxtargntarg)./(r2) +...
            1/(pi).*(rn.*rtautarg.*tauxtargntarg)./(r2.^2) + ...
            1/(pi).*(rtautarg.*tauxtargnsrc.*rntarg)./(r2.^2) -....
            2/(pi).*(rn.*rntarg.*(rtautarg.*rtautarg))./(r2.^3) + ...
            1/(2*pi).*(rn.*rntarg)./(r2.^2)) - 3/(4*pi).*normaldot./r2 +...
            ((2-coefs(1))/(2*pi)).*((rtautarg.*rtautarg).*(normaldot))./(r2.^2) -...
            ((2-coefs(1))/(4*pi)).*(normaldot)./r2 - ...
          ((1+coefs(1))/(4*pi)).*((tangentdot)./(r2) - 2*(rtautarg).*(rtau)./(r2.^2)) + temp + eulerconstantpart + puresmoothpart;
   end
end

% kernels in K11 with hilbert transform subtractions. 
% (i.e. beta*(G_{nx nx tauy} + 1/4 H + nu*(G_{taux taux tauy} + 1/4 H))

if strcmpi(type, 'free plate hilbert subtract')                                 
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};

   [~, third, ~] = chnk.helm2d.thirdforth_derivatives(zk, src, targ);            % Hankel part
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   zkimag = (1i)*zk;
   [~, thirdK, ~] = chnk.helm2d.thirdforth_derivatives(zkimag, src, targ);     % modified bessel K part
   [~,grad] = chnk.lap2d.green(src,targ,true); 

  
  

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;
    

   hilb = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx);

   submat = ((1+ coefs(1))/2)*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 4).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 6).*(nytarg.*nytarg.*tauy)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(nxtarg.*nxtarg.*taux) + thirdK(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        thirdK(:, :, 4).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        thirdK(:, :, 6).*(nytarg.*nytarg.*tauy)) + 0.25*hilb) +...
       ((1+ coefs(1))/2)*coefs(1).*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 4).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 6).*(tauytarg.*tauytarg.*tauy)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + thirdK(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       thirdK(:, :, 4).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        thirdK(:, :, 6).*(tauytarg.*tauytarg.*tauy)) + 0.25*hilb);   % hilbert subtraction 

end

% kernels in K21 that are coupled with Hilbert transforms. 
% (i.e beta* G_{nx nx nx tauy} + (2-nu)*beta*G_{taux taux tauy  nx})
if strcmpi(type, 'free plate coupled hilbert')                                  
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};

   [~, ~, forth] = chnk.helm2d.thirdforth_derivatives(zk, src, targ);            % Hankel part
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   zkimag = (1i)*zk;
   [~, ~, forthK] = chnk.helm2d.thirdforth_derivatives(zkimag, src, targ);     % modified bessel K part
  

  
  

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;
    


   submat =  ((1+ coefs(1))/2).*(1/(2*zk^2).*(forth(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*taux) + ...
          forth(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*tauy + 3*nxtarg.*nxtarg.*nytarg.*taux) + ...
          forth(:, :, 4).*(3*nxtarg.*nxtarg.*nytarg.*tauy + 3*nxtarg.*nytarg.*nytarg.*taux) +...
          forth(:, :, 6).*(3*nxtarg.*nytarg.*nytarg.*tauy + nytarg.*nytarg.*nytarg.*taux) +...
          forth(:, :, 8).*(nytarg.*nytarg.*nytarg.*tauy)) - ...
          1/(2*zk^2).*(forthK(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*taux) + ...
          forthK(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*tauy + 3*nxtarg.*nxtarg.*nytarg.*taux) + ...
          forthK(:, :, 4).*(3*nxtarg.*nxtarg.*nytarg.*tauy + 3*nxtarg.*nytarg.*nytarg.*taux) +...
          forthK(:, :, 6).*(3*nxtarg.*nytarg.*nytarg.*tauy + nytarg.*nytarg.*nytarg.*taux) +...
          forthK(:, :, 8).*(nytarg.*nytarg.*nytarg.*tauy))) +...
          ((2-coefs(1))/2)*(1+coefs(1)).*(1/(2*zk^2).*(forth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*taux) + ...
          forth(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*tauy + tauxtarg.*tauxtarg.*nytarg.*taux + 2*tauxtarg.*tauytarg.*nxtarg.*taux) + ...
          forth(:, :, 4).*(tauxtarg.*tauxtarg.*nytarg.*tauy + 2*tauxtarg.*tauytarg.*nxtarg.*tauy + tauytarg.*tauytarg.*nxtarg.*taux + 2*tauxtarg.*tauytarg.*nytarg.*taux) +...
          forth(:, :, 6).*(tauytarg.*tauytarg.*nxtarg.*tauy + 2*tauxtarg.*tauytarg.*nytarg.*tauy + tauytarg.*tauytarg.*nytarg.*taux) +...
         forth(:, :, 8).*(tauytarg.*tauytarg.*nytarg.*tauy)) - ...
         1/(2*zk^2).*(forthK(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*taux) + ...
         forthK(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*tauy + tauxtarg.*tauxtarg.*nytarg.*taux + 2*tauxtarg.*tauytarg.*nxtarg.*taux) + ...
         forthK(:, :, 4).*(tauxtarg.*tauxtarg.*nytarg.*tauy + 2*tauxtarg.*tauytarg.*nxtarg.*tauy + tauytarg.*tauytarg.*nxtarg.*taux + 2*tauxtarg.*tauytarg.*nytarg.*taux) +...
         forthK(:, :, 6).*(tauytarg.*tauytarg.*nxtarg.*tauy + 2*tauxtarg.*tauytarg.*nytarg.*tauy + tauytarg.*tauytarg.*nytarg.*taux) +...
         forthK(:, :, 8).*(tauytarg.*tauytarg.*nytarg.*tauy)));
    
end

% Updated part in K21 that is not coupled with Hilbert transform. 
%(i.e. (1-nu)*(-G_{nx nx ny} + G_{taux taux ny}). 

if strcmpi(type, 'free plate K21 second part')                                                         
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};

    
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;
   
   zkimag = (1i)*zk;

   [~, third, ~] = chnk.helm2d.thirdforth_derivatives(zk, src, targ);            % Hankel part
   [~, thirdK, ~] = chnk.helm2d.thirdforth_derivatives(zkimag, src, targ);     % modified bessel K part

   submat = -((1-coefs(1))/(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nx) + third(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
        third(:, :, 4).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
        third(:, :, 6).*(nytarg.*nytarg.*ny)) - ...
        (1-coefs(1))/(2*zk^2).*(thirdK(:, :, 1).*(nxtarg.*nxtarg.*nx) + thirdK(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
        thirdK(:, :, 4).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
        thirdK(:, :, 6).*(nytarg.*nytarg.*ny))) + ...
       (1-coefs(1))/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + third(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
        third(:, :, 4).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
        third(:, :, 6).*(tauytarg.*tauytarg.*ny)) - ...
        (1-coefs(1))/(2*zk^2).*(thirdK(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + thirdK(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
        thirdK(:, :, 4).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
        thirdK(:, :, 6).*(tauytarg.*tauytarg.*ny));                                 % (1-nu)*(- G_{nx nx ny} + G_{ny taux taux})  
end

% Updated part in K21 that is coupled with Hilbert transform. 
%(i.e. (1-nu)*(beta*(G_{taux taux tauy} + 1/4 H) - beta*(G_{nx nx tauy} + 1/4 H)). 

if strcmpi(type, 'free plate K21 hilbert part')                                                         
  srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};

   [~, third, ~] = chnk.helm2d.thirdforth_derivatives(zk, src, targ);            % Hankel part
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   zkimag = (1i)*zk;
   [~, thirdK, ~] = chnk.helm2d.thirdforth_derivatives(zkimag, src, targ);     % modified bessel K part
   [~,grad] = chnk.lap2d.green(src,targ,true); 

  
  

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;
    

   hilb = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx);

   submat =  (1-coefs(1)).*(((1+ coefs(1))/2).*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 4).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 6).*(tauytarg.*tauytarg.*tauy)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + thirdK(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       thirdK(:, :, 4).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        thirdK(:, :, 6).*(tauytarg.*tauytarg.*tauy)) + 0.25*hilb) - ...
        ((1+ coefs(1))/2)*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 4).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 6).*(nytarg.*nytarg.*tauy)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(nxtarg.*nxtarg.*taux) + thirdK(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        thirdK(:, :, 4).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        thirdK(:, :, 6).*(nytarg.*nytarg.*tauy)) + 0.25*hilb)) ;



end

% Updated part in K22. (i.e. (1-nu)*(G_{taux taux} - G_{nx nx})

if strcmpi(type, 'free plate K22 second part')
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};



   [hess, ~, ~] = chnk.helm2d.thirdforth_derivatives(zk, src, targ);            % Hankel part
   
   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   zkimag = (1i)*zk;
   [hessK, ~, ~] = chnk.helm2d.thirdforth_derivatives(zkimag, src, targ);     % modified bessel K part
   
   

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 


   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;
    
   
   submat =  (1-coefs(1)).*(1/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg)-...
           1/(2*zk^2).*(hessK(:, :, 1).*tauxtarg.*tauxtarg + hessK(:, :, 2).*(2*tauxtarg.*tauytarg) + ...
           hessK(:, :, 3).*tauytarg.*tauytarg)-...
           (1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg)-...
           1/(2*zk^2).*(hessK(:, :, 1).*nxtarg.*nxtarg + hessK(:, :, 2).*(2*nxtarg.*nytarg) + hessK(:, :, 3).*nytarg.*nytarg)));


end







if strcmpi(type, 'test kernel')                                                         % test G_{taux taux tauy nx}
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};


   [hess, third, forth] = chnk.helm2d.thirdforth_derivatives(zk, src, targ);           % Hankel part
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   
   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   zkimag = (1i)*zk;
   [hessK, thirdK, forthK] = chnk.helm2d.thirdforth_derivatives(zkimag, src, targ);     % modified bessel K part
   [~,grad] = chnk.lap2d.green(src,targ,true);

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;


   hilb = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx);

   submat = 1/(2*zk^2).*(forth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*taux) + ...
          forth(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*tauy + tauxtarg.*tauxtarg.*nytarg.*taux + 2*tauxtarg.*tauytarg.*nxtarg.*taux) + ...
          forth(:, :, 4).*(tauxtarg.*tauxtarg.*nytarg.*tauy + 2*tauxtarg.*tauytarg.*nxtarg.*tauy + tauytarg.*tauytarg.*nxtarg.*taux + 2*tauxtarg.*tauytarg.*nytarg.*taux) +...
          forth(:, :, 6).*(tauytarg.*tauytarg.*nxtarg.*tauy + 2*tauxtarg.*tauytarg.*nytarg.*tauy + tauytarg.*tauytarg.*nytarg.*taux) +...
         forth(:, :, 8).*(tauytarg.*tauytarg.*nytarg.*tauy)) - ...
         1/(2*zk^2).*(forthK(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*taux) + ...
         forthK(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*tauy + tauxtarg.*tauxtarg.*nytarg.*taux + 2*tauxtarg.*tauytarg.*nxtarg.*taux) + ...
         forthK(:, :, 4).*(tauxtarg.*tauxtarg.*nytarg.*tauy + 2*tauxtarg.*tauytarg.*nxtarg.*tauy + tauytarg.*tauytarg.*nxtarg.*taux + 2*tauxtarg.*tauytarg.*nytarg.*taux) +...
         forthK(:, :, 6).*(tauytarg.*tauytarg.*nxtarg.*tauy + 2*tauxtarg.*tauytarg.*nytarg.*tauy + tauytarg.*tauytarg.*nytarg.*taux) +...
         forthK(:, :, 8).*(tauytarg.*tauytarg.*nytarg.*tauy));
end                                                         





% first-kernel of evaluation for the clamped plate
if strcmpi(type, 'first kernel')
    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    dx = repmat(srctang(1,:),nt,1);
    dy = repmat(srctang(2,:),nt,1);
    ds = sqrt(dx.*dx+dy.*dy);

    taux = dx./ds;
    tauy = dy./ds;



    [~, third, ~] = chnk.helm2d.thirdforth_derivatives(zk, src, targ);           % Hankel part
    
    zkimag = 1i*zk;
    [~, thirdK, ~] = chnk.helm2d.thirdforth_derivatives(zkimag, src, targ);     % modified bessel K part



    submat = 1/(2*zk^2).*(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
       third(:, :, 4).*(3*nx.*ny.*ny) + third(:, :, 6).*(ny.*ny.*ny)) - ...
        1/(2*zk^2).*(thirdK(:, :, 1).*(nx.*nx.*nx) + thirdK(:, :, 2).*(3*nx.*nx.*ny) +...
       thirdK(:, :, 4).*(3*nx.*ny.*ny) + thirdK(:, :, 6).*(ny.*ny.*ny)) + ...
       (3/(2*zk^2).*(third(:, :, 1).*(nx.*taux.*taux) + third(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       third(:, :, 4).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + third(:, :, 6).*(ny.*tauy.*tauy))-...
       3/(2*zk^2).*(thirdK(:, :, 1).*(nx.*taux.*taux) + thirdK(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       thirdK(:, :, 4).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + thirdK(:, :, 6).*(ny.*tauy.*tauy)));

end

% second-kernel of evaluation for the clamped plate
if strcmpi(type, 'second kernel')
    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    dx = repmat(srctang(1,:),nt,1);
    dy = repmat(srctang(2,:),nt,1);

    ds = sqrt(dx.*dx+dy.*dy);

    taux = dx./ds;
    tauy = dy./ds;

    [hess, ~, ~] = chnk.helm2d.thirdforth_derivatives(zk, src, targ);           % Hankel part
    
    zkimag = 1i*zk;
    [hessK, ~, ~] = chnk.helm2d.thirdforth_derivatives(zkimag, src, targ);     % modified bessel K part

    submat =  -(1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny))-...
           1/(2*zk^2).*(hessK(:, :, 1).*(nx.*nx) + hessK(:, :, 2).*(2*nx.*ny) + hessK(:, :, 3).*(ny.*ny)))+...
          (1/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy))-...
           1/(2*zk^2).*(hessK(:, :, 1).*(taux.*taux) + hessK(:, :, 2).*(2*taux.*tauy) + hessK(:, :, 3).*(tauy.*tauy))); % -G_{ny ny}  + G_{tauy tauy}
end


if strcmpi(type, 'free plate eval first')                                               % G_{ny}
   srcnorm = srcinfo.n;
   [~,grad] = chnk.helm2d.green(zk,src,targ);        % Hankel part
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   


   zkimag = (1i)*zk;
   [~,gradK] = chnk.helm2d.green(zkimag,src,targ);    % modified bessel K part

 
  

   submat = (-1/(2*zk^2).*(grad(:, :, 1).*(nx) + grad(:, :, 2).*ny) + ...
       1/(2*zk^2).*(gradK(:, :, 1).*(nx) + gradK(:, :, 2).*ny)); 


end

if strcmpi(type, 'free plate eval first hilbert')                                               % G_{tauy}
   srctang = srcinfo.d;
   coefs = varargin{1};

   [~,grad] = chnk.helm2d.green(zk,src,targ);        % Hankel part

  
   zkimag = (1i)*zk;
   [~,gradK] = chnk.helm2d.green(zkimag,src,targ);    % modified bessel K part

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   ds = sqrt(dx.*dx+dy.*dy);
 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;



    submat = ((1 + coefs(1))/2).*(-1/(2*zk^2).*(grad(:, :, 1).*(taux) + grad(:, :, 2).*tauy) + ...
       1/(2*zk^2).*(gradK(:, :, 1).*(taux) + gradK(:, :, 2).*tauy));                    % G_{tauy}
end


if strcmpi(type, 'free plate eval second')                                          % G = 1/(2k^2) (i/4 H_0^{1} - 1/2pi K_0)

   [val,~] = chnk.helm2d.green(zk,src,targ);        % Hankel part

   zkimag = (1i)*zk;
   [valK,~] = chnk.helm2d.green(zkimag,src,targ);    % modified bessel K part

   submat = 1/(2*zk^2).*val - 1/(2*zk^2).*valK;

end






% Dirichlet data/potential correpsonding to transmission rep
if strcmpi(type,'trans_rep') 

  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end;
  srcnorm = srcinfo.n;
  [submats,grad] = chnk.helm2d.green(zk,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  submatd = -(grad(:,:,1).*nx + grad(:,:,2).*ny);

  submat = zeros(1*nt,2*ns);
  submat(1:1:1*nt,1:2:2*ns) = coef(1)*submatd;
  submat(1:1:1*nt,2:2:2*ns) = coef(2)*submats;
end

% Neumann data corresponding to transmission rep
if strcmpi(type,'trans_rep_prime')
  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end;
  targnorm = targinfo.n;
  srcnorm = srcinfo.n;
  
  submat = zeros(nt,ns);

  % Get gradient and hessian info
  [submats,grad,hess] = chnk.helm2d.green(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);

  % D'
  submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
      + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
  % S'
  submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
  submat = zeros(nt,2*ns)
  submat(1:1:1*nt,1:2:2*ns) = coef(1)*submatdp;
  submat(1:1:1*nt,2:2:2*ns) = coef(2)*submatsp;
end


% Gradient correpsonding to transmission rep
if strcmpi(type,'trans_rep_grad')
  coef = ones(2,1);
  if(nargin == 5); coef = varargin{1}; end;
  
  srcnorm = srcinfo.n;
  
  submat = zeros(nt,ns,6);
  % S
  [submats,grad,hess] = chnk.helm2d.green(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  % D
  submatd  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);

  submat = zeros(2*nt,2*ns);
  
  submat(1:2:2*nt,1:2:2*ns) = -coef(1)*(hess(:,:,1).*nxsrc + hess(:,:,2).*nysrc);
  submat(1:2:2*nt,2:2:2*ns) = coef(2)*grad(:,:,1);
    
  submat(2:2:2*nt,1:2:2*ns) = -coef(1)*(hess(:,:,2).*nxsrc + hess(:,:,3).*nysrc);
  submat(2:2:2*nt,2:2:2*ns) = coef(2)*grad(:,:,2);








end


