function submat= kern(zk,srcinfo,targinfo,type,varargin)
%FLEX2D.KERN flexural wave kernels in 2D
% 
% Syntax: submat = chnk.flex2d.kern(zk,srcinfo,targingo,type,varargin)
%
% Let x be targets and y be sources for these formulas, with
% n_x and n_y the corresponding unit normals at those points
% (if defined). Note that the normal information is obtained
% by taking the perpendicular to the provided tangential deriviative
% info and normalizing  
%  
% Kernels based on the flexural wave Green's function:
%         G(x,y)=1/zk^2*(i/4 H_0^(1)(k|x-y|) - 1/(2 pi)*K_0(k|x-y|))
%
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
%                be provided. sprime requires normal info.
%   type - string, determines kernel type
%                type == 's', single layer kernel S
%                type == 'sprime', normal derivative of single
%                      layer S'
%                type == 'clamped_plate_bcs', returns the clamped plate
%                         boundary conditions applied to a point source.
%                         requires target normal info.
%                type == 'clamped_plate', returns the 2x2 system of kernels
%                         needed to solve the clamped plate problem, 
%                          [G_{n_y n_y n_y}+3*G_{n_y tau_y tau_y}, 
%                           -G_{n_y n_y}+G_{tau_y tau_y}; 
%                            G_{n_x n_y n_y n_y}+3*G_{n_x n_y tau_y tau_y}, 
%                           -G_{n_x n_y n_y}+G_{n_x tau_y tau_y}]
%                type == 'clamped_plate_eval', returns the kernels
%                        corresponding to the clamped plate representation, 
%                        [G_{n_y n_y n_y}+3*G_{n_y tau_y tau_y}, 
%                         -G_{n_y n_y}+G_{tau_y tau_y}]
%                type == 'free_plate_bcs', returns the free plate
%                         boundary conditions applied to a point source.
%                         requires target normals, tangents, and d2.
%                type == 'free_plate', returns the 4x2 system of kernels 
%                         needed to solve the free plate integral equation,
%                         [K11, K12; K21, K22; K11H, 0; K21H, 0] where 
%                         K11H and K21H are kernels that need to be coupled
%                         with the Hilbert transform. For more info see the
%                         reference below. 
%                type == 'free_plate_eval' returns the kernels
%                         corresponding to the free plate representation:
%                         [G_{n_y}, beta*G_{tau_y}, G]
%                type == 'supported_plate_bcs', returns the supported plate
%                         boundary conditions applied to a point source.
%                         requires target normal and tangent info
%                type == 'supported_plate_log', returns the 2x2 system of 
%                         kernels for the supported plate integral 
%                         equation, [K11, K12; K21, K22] that have to be 
%                         discretized using log quads. this kernel requires 
%                         kappa' in the first data field of srcinfo. For 
%                         more info see the reference below. 
%                type == 'supported_plate_smooth', returns the one kernel 
%                         in the supported plate equation that has to be 
%                         discretized using smooth quads to avoid close 
%                         evaluations. this kernel requires kappa' in the 
%                         first data field of srcinfo and kappa'' in the
%                         second data field.
%                type == 'supported_plate_eval' returns the kernels
%                         corresponding to the supported plate
%                         representation. also requires 
%                         kappa' in the first data field of srcinfo.
%
%   varargin{1} - nu: scalar value for the Poisson's ratio of the plate.
%                    Needed whenever solving the supported and free plates, 
%                    not needed for the clamped plate
%
% Output:
%   submat - the evaluation of the selected kernel(s) for the
%            provided sources and targets
%
% For more info on the integral equation methods implemented here see:
% Nekrasov, P., Su, Z., Askham, T., & Hoskins, J. G. (2024). Boundary 
% Integral Formulations for Flexural Wave Scattering in Thin Plates. arXiv 
% preprint arXiv:2409.19160.
% 
% see also CHNK.FLEX2D.HKDIFFGREEN
  
src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);

%%% STANDARD LAYER POTENTIALS

switch lower(type)
case {'s', 'single'} % flexural wave single layer

   val = chnk.flex2d.hkdiffgreen(zk,src,targ);  
   submat = 1/(2*zk^2).*val;

case {'sp', 'sprime'} % normal derivative of flexural wave single layer

   targnorm = targinfo.n;
   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   [~,grad] = chnk.flex2d.hkdiffgreen(zk,src,targ);  
   submat = 1/(2*zk^2).*(grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);

%%% CLAMPED PLATE KERNELS

% boundary conditions applied to a point source
case {'clamped_plate_bcs'}
    nxtarg = targinfo.n(1,:).'; 
    nytarg = targinfo.n(2,:).';  
    submat = zeros(2*nt,ns);
    
    [val, grad] = chnk.flex2d.hkdiffgreen(zk, src, targ);
    
    firstbc = 1/(2*zk^2).*val ;
    secondbc = 1/(2*zk^2).*(grad(:, :, 1).*nxtarg + grad(:, :, 2).*nytarg);
   
    submat(1:2:end,:) = firstbc;
    submat(2:2:end,:) = secondbc;

% kernels for the clamped plate integral equation
case {'clamped_plate'}
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   
   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);
   
   [~, ~, hess, third, ~] = chnk.flex2d.hkdiffgreen(zk, src, targ); 
   [~, ~, ~, ~, fourth] = chnk.flex2d.hkdiffgreen(zk, src, targ, true);

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);
    
   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds;
   tauy = dy./ds;
   
   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;

   rn = rx.*nx + ry.*ny;
   rtau = rx.*taux + ry.*tauy;
   ntargtau = nxtarg.*taux + nytarg.*tauy;

   rntarg = rx.*nxtarg + ry.*nytarg;

   K11 = -(1/(2*zk^2).*(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
       third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny))) - ...
       (3/(2*zk^2).*(third(:, :, 1).*(nx.*taux.*taux) + third(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       third(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + third(:, :, 4).*(ny.*tauy.*tauy)));  % G_{ny ny ny} + 3G_{ny tauy tauy}

   K12 = -(1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny)))+...
          (1/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy))); % -G_{ny ny}  + G_{tauy tauy}

   K21 = -(1/(2*zk^2).*(fourth(:, :, 1).*(nx.*nx.*nx.*nxtarg) + fourth(:, :, 2).*(nx.*nx.*nx.*nytarg + 3*nx.*nx.*ny.*nxtarg) + ...
          fourth(:, :, 3).*(3*nx.*nx.*ny.*nytarg + 3*nx.*ny.*ny.*nxtarg) + fourth(:, :, 4).*(3*nx.*ny.*ny.*nytarg +ny.*ny.*ny.*nxtarg)+...
          fourth(:, :, 5).*(ny.*ny.*ny.*nytarg)) ) - ...
          (3/(2*zk^2).*(fourth(:, :, 1).*(nx.*taux.*taux.*nxtarg)+ fourth(:, :, 2).*(nx.*taux.*taux.*nytarg + 2*nx.*taux.*tauy.*nxtarg + ny.*taux.*taux.*nxtarg) +...
          fourth(:, :, 3).*(2*nx.*taux.*tauy.*nytarg + ny.*taux.*taux.*nytarg + nx.*tauy.*tauy.*nxtarg + 2*ny.*taux.*tauy.*nxtarg) + ...
          fourth(:, :, 4).*(nx.*tauy.*tauy.*nytarg +2*ny.*taux.*tauy.*nytarg + ny.*tauy.*tauy.*nxtarg) +...
          fourth(:, :, 5).*(ny.*tauy.*tauy.*nytarg))) + ...
          1/pi.*(-3*rn.*rntarg./(r2.^2) + 4.*(rn.^3).*rntarg./(r2.^3) + 3*(rn.*rtau.*ntargtau)./ (r2.^2));

   K22 = -(1/(2*zk^2).*(third(:,:, 1).*(nx.*nx.*nxtarg) +third(:, :, 2).*(nx.*nx.*nytarg + 2*nx.*ny.*nxtarg) + third(:, :, 3).*(2*nx.*ny.*nytarg + ny.*ny.*nxtarg)+...
         third(:, :,4).*(ny.*ny.*nytarg))) + ...
         (1/(2*zk^2).*(third(:,:, 1).*(taux.*taux.*nxtarg) +third(:, :, 2).*(taux.*taux.*nytarg + 2*taux.*tauy.*nxtarg) + third(:, :, 3).*(2*taux.*tauy.*nytarg + tauy.*tauy.*nxtarg)+...
         third(:, :,4).*(tauy.*tauy.*nytarg)));

  submat = zeros(2*nt,2*ns);
  
  submat(1:2:end,1:2:end) = K11;
  submat(1:2:end,2:2:end) = K12;
    
  submat(2:2:end,1:2:end) = K21;
  submat(2:2:end,2:2:end) = K22;

% clamped plate kernels for plotting
case {'clamped_plate_eval'}

    submat = zeros(nt,2*ns);

    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    dx = repmat(srctang(1,:),nt,1);
    dy = repmat(srctang(2,:),nt,1);
    ds = sqrt(dx.*dx+dy.*dy);

    taux = dx./ds;
    tauy = dy./ds;

    [~, ~, hess, third] = chnk.flex2d.hkdiffgreen(zk, src, targ);           % Hankel part

    K1 = -(1/(2*zk^2).*(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
       third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny)) ) - ...
       (3/(2*zk^2).*(third(:, :, 1).*(nx.*taux.*taux) + third(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       third(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + third(:, :, 4).*(ny.*tauy.*tauy)));  % G_{ny ny ny} + 3G_{ny tauy tauy}

    K2 =  -(1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny)))+...
          (1/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy))); % -G_{ny ny}  + G_{tauy tauy}

    submat(:,1:2:end) = K1;
    submat(:,2:2:end) = K2;


%%% FREE PLATE KERNELS

% boundary conditions applied to a point source
case {'free_plate_bcs'}
    targnorm = targinfo.n;
    targtang = targinfo.d;
    targd2 = targinfo.d2;
    nu = varargin{1};
    
    [~, ~, hess, third] = chnk.flex2d.hkdiffgreen(zk, src, targ);

    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
    
    dx1 = repmat((targtang(1,:)).',1,ns);
    dy1 = repmat((targtang(2,:)).',1,ns);
    
    ds1 = sqrt(dx1.*dx1+dy1.*dy1); 
    
    d2x1 = repmat((targd2(1,:)).',1,ns);
    d2y1 = repmat((targd2(2,:)).',1,ns);
    
    tauxtarg = dx1./ds1;
    tauytarg = dy1./ds1;
    
    denom = sqrt(dx1.^2 + dy1.^2).^3;
    numer = dx1.*d2y1 - d2x1.*dy1;
    
    kappatarg = numer ./ denom; % target curvature
    
    firstbc = 1/(2*zk^2).*(hess(:, :, 1).*(nxtarg.*nxtarg) + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*(nytarg.*nytarg))+...
    nu/(2*zk^2).*(hess(:, :, 1).*(tauxtarg.*tauxtarg) + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*(tauytarg.*tauytarg));
    
    secondbc = 1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nxtarg) + third(:, :, 2).*(3*nxtarg.*nxtarg.*nytarg) +...
    third(:, :, 3).*(3*nxtarg.*nytarg.*nytarg) + third(:, :, 4).*(nytarg.*nytarg.*nytarg))+...
    (2-nu)/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg) + third(:, :, 2).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
    third(:, :, 3).*(2*tauxtarg.*tauytarg.*nytarg+ tauytarg.*tauytarg.*nxtarg) +...
    + third(:, :, 4).*(tauytarg.*tauytarg.*nytarg))+...
    (1-nu).*kappatarg.*(1/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg)-...
    (1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg)));

    submat = zeros(2*nt,ns);
    submat(1:2:end,:) = firstbc;
    submat(2:2:end,:) = secondbc;

% kernels for the free plate integral equation 
case {'free_plate'}
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   targd2 = targinfo.d2;
   nu = varargin{1};

   [~, ~, hess, third, fourth] = chnk.flex2d.hkdiffgreen(zk, src, targ);     

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);


   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   ds = sqrt(dx.*dx+dy.*dy); 

   taux = dx ./ ds;
   tauy = dy ./ ds;

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);

   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   d2x1 = repmat((targd2(1,:)).',1,ns);
   d2y1 = repmat((targd2(2,:)).',1,ns);

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   denom = sqrt(dx1.^2 + dy1.^2).^3;
   numer = dx1.*d2y1 - d2x1.*dy1;

   kappatarg = numer ./ denom; % target curvature

   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;

   [~,grad] = chnk.lap2d.green(src,targ,true); 
   hilb = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx);
   
   K11 = -(1/(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nx) + third(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
        third(:, :, 4).*(nytarg.*nytarg.*ny))) - ...
       nu./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + third(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*ny)) ;  % first kernel with no hilbert transforms (G_{nx nx ny + nu G_{taux taux ny}).

   K12 =  1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg)+...
           nu/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg) ;    % G_{nx nx} + nu G_{taux taux}
   
   K21 = kappatarg./(2*zk^2).*(1-nu).*((third(:, :, 1).*(nxtarg.*nxtarg.*nx) + third(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
            third(:, :, 3).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
            third(:, :, 4).*(nytarg.*nytarg.*ny)) - ...
           (third(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + third(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
            third(:, :, 3).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
            third(:, :, 4).*(tauytarg.*tauytarg.*ny)) ) ...
        - 1/(2*zk^2).*(fourth(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*nx) + fourth(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*ny + 3*nxtarg.*nxtarg.*nytarg.*nx) + ...
          fourth(:, :, 3).*(3*nxtarg.*nxtarg.*nytarg.*ny + 3*nxtarg.*nytarg.*nytarg.*nx) +...
          fourth(:, :, 4).*(3*nxtarg.*nytarg.*nytarg.*ny +nytarg.*nytarg.*nytarg.*nx)+...
          fourth(:, :, 5).*(nytarg.*nytarg.*nytarg.*ny)) - ...
          ((2-nu)/(2*zk^2).*(fourth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*nx) + fourth(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*ny + tauxtarg.*tauxtarg.*nytarg.*nx + 2*tauxtarg.*tauytarg.*nxtarg.*nx)+...
          fourth(:, :, 3).*(tauxtarg.*tauxtarg.*nytarg.*ny + 2*tauxtarg.*tauytarg.*nxtarg.*ny + tauytarg.*tauytarg.*nxtarg.*nx + 2*tauxtarg.*tauytarg.*nytarg.*nx)+...
          fourth(:, :, 4).*(tauytarg.*tauytarg.*nxtarg.*ny + 2*tauxtarg.*tauytarg.*nytarg.*ny + tauytarg.*tauytarg.*nytarg.*nx) +...
          fourth(:, :, 5).*(tauytarg.*tauytarg.*nytarg.*ny)) ) - ...          
          (1+nu)/(4*pi).*((taux.*tauxtarg + tauy.*tauytarg)./(r2) - 2*(rx.*tauxtarg + ry.*tauytarg).*(rx.*taux + ry.*tauy)./(r2.^2));

   K22 = 1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nxtarg) + third(:, :, 2).*(3*nxtarg.*nxtarg.*nytarg) +...
       third(:, :, 3).*(3*nxtarg.*nytarg.*nytarg) + third(:, :, 4).*(nytarg.*nytarg.*nytarg))  +...
        (2-nu)/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg) + third(:, :, 2).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*nytarg + tauytarg.*tauytarg.*nxtarg) +...
        + third(:, :, 4).*(tauytarg.*tauytarg.*nytarg)) + ... % G_{nx nx nx} + (2-nu) G_{taux taux nx}
        + kappatarg.*(1-nu).*(1/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg)-...
           1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg));
    
   K11H =  -(1+ nu)/2*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 4).*(nytarg.*nytarg.*tauy)) )  - ...
       (1+ nu)/2*nu.*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy))) + (1+ nu)/2*(1+nu).*0.25*hilb ;

   K21H = kappatarg.*(1-nu).*(-((1+ nu)/2).*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy)))  + ...
        ((1+ nu)/2)*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 4).*(nytarg.*nytarg.*tauy)))) ...
        -(1+ nu)/2.*(1/(2*zk^2).*(fourth(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*taux) + ...
          fourth(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*tauy + 3*nxtarg.*nxtarg.*nytarg.*taux) + ...
          fourth(:, :, 3).*(3*nxtarg.*nxtarg.*nytarg.*tauy + 3*nxtarg.*nytarg.*nytarg.*taux) +...
          fourth(:, :, 4).*(3*nxtarg.*nytarg.*nytarg.*tauy + nytarg.*nytarg.*nytarg.*taux) +...
          fourth(:, :, 5).*(nytarg.*nytarg.*nytarg.*tauy)) ) - ...
          ((2-nu)/2)*(1+nu).*(1/(2*zk^2).*(fourth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*taux) + ...
          fourth(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*tauy + tauxtarg.*tauxtarg.*nytarg.*taux + 2*tauxtarg.*tauytarg.*nxtarg.*taux) + ...
          fourth(:, :, 3).*(tauxtarg.*tauxtarg.*nytarg.*tauy + 2*tauxtarg.*tauytarg.*nxtarg.*tauy + tauytarg.*tauytarg.*nxtarg.*taux + 2*tauxtarg.*tauytarg.*nytarg.*taux) +...
          fourth(:, :, 4).*(tauytarg.*tauytarg.*nxtarg.*tauy + 2*tauxtarg.*tauytarg.*nytarg.*tauy + tauytarg.*tauytarg.*nytarg.*taux) +...
         fourth(:, :, 5).*(tauytarg.*tauytarg.*nytarg.*tauy)) ) ;
    
    submat = zeros(4*nt,2*ns);
    
    submat(1:4:end,1:2:end) = K11;
    submat(1:4:end,2:2:end) = K12;
    
    submat(2:4:end,1:2:end) = K21;
    submat(2:4:end,2:2:end) = K22;
    
    submat(3:4:end,1:2:end) = K11H;
    submat(4:4:end,1:2:end) = K21H;

% free plate kernels used for plotting 
case {'free_plate_eval'}
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   nu = varargin{1};

   [val,grad] = chnk.flex2d.hkdiffgreen(zk,src,targ); 
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds; 
   tauy = dy./ds;
  
   K1 = (-1/(2*zk^2).*(grad(:, :, 1).*(nx) + grad(:, :, 2).*ny)); 
   K1H = ((1 + nu)/2).*(-1/(2*zk^2).*(grad(:, :, 1).*(taux) + grad(:, :, 2).*tauy));                    % G_{tauy}
   K2 = 1/(2*zk^2).*val;

   submat = zeros(nt,3*ns);
   submat(:,1:3:end) = K1;
   submat(:,2:3:end) = K1H;
   submat(:,3:3:end) = K2;

%%% SUPPORTED PLATE KERNELS

% boundary conditions applied to a point source
case {'supported_plate_bcs'}
    nxtarg = targinfo.n(1,:).'; 
    nytarg = targinfo.n(2,:).';  
    dx = targinfo.d(1,:).';
    dy = targinfo.d(2,:).';
    ds = sqrt(dx.*dx+dy.*dy);
    tauxtarg = (dx./ds);                                                                       % normalization
    tauytarg = (dy./ds);

    nu = varargin{1};
    
    [val, ~, hess] = chnk.flex2d.hkdiffgreen(zk, src, targ);
    
    firstbc = 1/(2*zk^2).*val ;
    
    secondbc = 1/(2*zk^2).*(hess(:, :, 1).*(nxtarg.*nxtarg) + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*(nytarg.*nytarg))+...
               nu/(2*zk^2).*(hess(:, :, 1).*(tauxtarg.*tauxtarg) + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*(tauytarg.*tauytarg));
    
    submat = zeros(2*nt,ns);

    submat(1:2:end,:) = firstbc;
    submat(2:2:end,:) = secondbc;

% kernels for the supported plate integral equation that have to be 
% discretized using log quadrature
case {'supported_plate_log'}
    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    srcd2 = srcinfo.d2;
    targnorm = targinfo.n;
    targtang = targinfo.d;
    coefs = varargin{1};
    nu = coefs(1);
    
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    
    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
        
    dx = repmat(srctang(1,:),nt,1);
    dy = repmat(srctang(2,:),nt,1);
    
    ds = sqrt(dx.*dx+dy.*dy);
    
    taux = dx./ds;
    tauy = dy./ds;
    
    dx1 = repmat((targtang(1,:)).',1,ns);
    dy1 = repmat((targtang(2,:)).',1,ns);
    
    ds1 = sqrt(dx1.*dx1+dy1.*dy1); 
    
    tauxtarg = dx1./ds1;
    tauytarg = dy1./ds1;
    
    dx1 = repmat(srctang(1,:),nt,1);
    dy1 = repmat(srctang(2,:),nt,1);
    
    d2x1 = repmat(srcd2(1,:),nt,1);
    d2y1 = repmat(srcd2(2,:),nt,1);
    
    denom = sqrt(dx1.^2+dy1.^2).^3;
    numer = dx1.*d2y1-d2x1.*dy1;
    
    kappa = numer./denom;   
    
    kp = repmat(srcinfo.data(1,:),nt,1);
    
    a1 = 2-nu;
    a2 = (-1+nu)*(7+nu)/(3 - nu);
    a3 = (1-nu)*(3+nu)/(1+nu);
    
    [~, grad, hess, third, ~, ~] = chnk.flex2d.hkdiffgreen(zk, src, targ, false);      
    
    K11 = -1/(2*zk^2)*(third(:,:,1).*nx.^3 + 3*third(:,:,2).*nx.^2.*ny + 3*third(:,:,3).*nx.*ny.^2 + third(:,:,4).*ny.^3) + ...
         -a1/(2*zk^2)*(third(:,:,1).*nx.*taux.^2 + third(:,:,2).*(ny.*taux.^2 + 2*nx.*taux.*tauy) + third(:,:,3).*(nx.*tauy.^2 + 2*ny.*taux.*tauy) + third(:,:,4).*ny.*tauy.^2) + ...
         a2*kappa./(2*zk^2).*(hess(:,:,1).*nx.^2 + 2*hess(:,:,2).*nx.*ny + hess(:,:,3).*ny.^2) + ...
         -a3*kp./(2*zk^2).*(grad(:,:,1).*taux + grad(:,:,2).*tauy);
    
    K12 = -1/(2*zk^2).*(grad(:, :, 1).*nx + grad(:, :, 2).*ny);
           
    K22 = -1/(2*zk^2).*(third(:,:,1).*(nx.*nxtarg.^2) + third(:,:,2).*(2*nx.*nxtarg.*nytarg + ny.*nxtarg.^2) + ...
         + third(:,:,3).*(nx.*nytarg.^2 + 2*ny.*nxtarg.*nytarg) + third(:,:,4).*(ny.*nytarg.^2))+...
         -nu./(2*zk^2).*(third(:,:,1).*(nx.*tauxtarg.^2) + third(:,:,2).*(2*nx.*tauxtarg.*tauytarg + ny.*tauxtarg.^2) + ...
         + third(:,:,3).*(nx.*tauytarg.^2 + 2*ny.*tauxtarg.*tauytarg) + third(:,:,4).*(ny.*tauytarg.^2)); 
    
    
    [~, ~, ~, third, forth, fifth] = chnk.flex2d.hkdiffgreen(zk, src, targ, true);      
    
    K21 = -1./(2*zk^2).*(fifth(:,:,1).*nxtarg.^2.*nx.^3 + ... % K21 with the biharmonic Green's function subtracted off (biharmonic part handled in a separate function)
          fifth(:,:,2).*(2*nxtarg.*nytarg.*nx.^3+3*nxtarg.^2.*nx.^2.*ny) + ...
          fifth(:,:,3).*(nytarg.^2.*nx.^3 + 6*nxtarg.*nytarg.*nx.^2.*ny+3.*nxtarg.^2.*nx.*ny.^2) + ...
          fifth(:,:,4).*(3*nytarg.^2.*nx.^2.*ny+6.*nxtarg.*nytarg.*nx.*ny.^2+nxtarg.^2.*ny.^3)+ ...
          fifth(:,:,5).*(nytarg.*ny.^2.*(3*nytarg.*nx+2*nxtarg.*ny))+ ...
          fifth(:,:,6).*nytarg.^2.*ny.^3) + ... % G_{nx nx ny ny ny}
         -nu./(2*zk^2).*(fifth(:,:,1).*tauxtarg.^2.*nx.^3 + ...
          fifth(:,:,2).*(2*tauxtarg.*tauytarg.*nx.^3+3*tauxtarg.^2.*nx.^2.*ny) + ...
          fifth(:,:,3).*(tauytarg.^2.*nx.^3 + 6*tauxtarg.*tauytarg.*nx.^2.*ny+3.*tauxtarg.^2.*nx.*ny.^2) + ...
          fifth(:,:,4).*(3*tauytarg.^2.*nx.^2.*ny+6.*tauxtarg.*tauytarg.*nx.*ny.^2+tauxtarg.^2.*ny.^3)+ ...
          fifth(:,:,5).*(tauytarg.*ny.^2.*(3*tauytarg.*nx+2*tauxtarg.*ny))+ ...
          fifth(:,:,6).*tauytarg.^2.*ny.^3) + ... % nu*G_{taux taux ny ny ny}
        -a1./(2*zk^2).*(fifth(:,:,1).*nxtarg.^2.*nx.*taux.^2 + ...
          fifth(:,:,2).*(2.*nx.*nxtarg.^2.*taux.*tauy + 2.*nx.*nxtarg.*nytarg.*taux.^2 + nxtarg.^2.*ny.*taux.^2) + ...
          fifth(:,:,3).*(nx.*nxtarg.^2.*tauy.^2 + 4*nx.*nxtarg.*nytarg.*taux.*tauy + 2*nxtarg.^2.*ny.*taux.*tauy + nx.*nytarg.^2.*taux.^2 + 2*nxtarg.*ny.*nytarg.*taux.^2) + ...
          fifth(:,:,4).*(2*nx.*nxtarg.*nytarg.*tauy.^2 + nxtarg.^2.*ny.*tauy.^2 + 2*nx.*nytarg.^2.*taux.*tauy + 4*nxtarg.*ny.*nytarg.*taux.*tauy + ny.*nytarg.^2.*taux.^2) + ...
          fifth(:,:,5).*(nytarg.*tauy.*(2*ny.*nytarg.*taux + 2*nxtarg.*ny.*tauy + nx.*nytarg.*tauy)) + ...
          fifth(:,:,6).*ny.*nytarg.^2.*tauy.^2) + ...  % G_{nx nx ny tauy tauy}
        -nu*a1./(2*zk^2).*(fifth(:,:,1).*tauxtarg.^2.*nx.*taux.^2 + ...
          fifth(:,:,2).*(ny.*taux.^2.*tauxtarg.^2 + 2.*nx.*taux.*tauxtarg.^2.*tauy +2.*nx.*taux.^2.*tauxtarg.*tauytarg) + ...
          fifth(:,:,3).*(nx.*tauxtarg.^2.*tauy.^2 + 4*nx.*tauxtarg.*tauytarg.*taux.*tauy + 2*tauxtarg.^2.*ny.*taux.*tauy + nx.*tauytarg.^2.*taux.^2 + 2*tauxtarg.*ny.*tauytarg.*taux.^2) + ...
          fifth(:,:,4).*(2*nx.*tauxtarg.*tauytarg.*tauy.^2 + tauxtarg.^2.*ny.*tauy.^2 + 2*nx.*tauytarg.^2.*taux.*tauy + 4*tauxtarg.*ny.*tauytarg.*taux.*tauy + ny.*tauytarg.^2.*taux.^2) + ...
          fifth(:,:,5).*(tauytarg.*tauy.*(2*ny.*tauytarg.*taux + 2*tauxtarg.*ny.*tauy + nx.*tauytarg.*tauy)) + ...
          fifth(:,:,6).*ny.*tauytarg.^2.*tauy.^2) + ... % nu*G_{taux taux ny tauy tauy}
         a2.*kappa./(2*zk^2).*(forth(:, :, 1).*(nxtarg.*nxtarg.*nx.*nx) + ...
          forth(:, :, 2).*(2*nx.^2.*nxtarg.*nytarg + 2*nx.*nxtarg.^2.*ny) + ...
          forth(:, :, 3).*(nxtarg.^2.*ny.^2 + 4*nx.*nxtarg.*ny.*nytarg + nx.^2.*nytarg.^2) +...
          forth(:, :, 4).*(2*ny.*nytarg.*(nxtarg.*ny + nx.*nytarg))+...
          forth(:, :, 5).*(nytarg.*nytarg.*ny.*ny)) + ... % kappa*G_{nx nx ny ny}
         nu*a2.*kappa./(2*zk^2).*(forth(:, :, 1).*(tauxtarg.*tauxtarg.*nx.*nx) + ...
          forth(:, :, 2).*(2*nx.^2.*tauxtarg.*tauytarg + 2*nx.*tauxtarg.^2.*ny) + ...
          forth(:, :, 3).*(tauxtarg.^2.*ny.^2 + 4*nx.*tauxtarg.*ny.*tauytarg + nx.^2.*tauytarg.^2) +...
          forth(:, :, 4).*(2*ny.*tauytarg.*(tauxtarg.*ny + nx.*tauytarg))+...
          forth(:, :, 5).*(tauytarg.*tauytarg.*ny.*ny)) + ... % kappa*nu*G_{nx nx ny ny}
        -a3.*kp./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
            third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
            third(:, :, 4).*(nytarg.*nytarg.*tauy)) + ... % kp*G_{nx nx tauy}
         -nu*a3.*kp./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
            third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
            third(:, :, 4).*(tauytarg.*tauytarg.*tauy)) ;  % kp*nu*G_{taux taux tauy}
    
    submat = zeros(2*nt,2*ns);
    
    submat(1:2:end,1:2:end) = K11;
    submat(1:2:end,2:2:end) = K12;
    
    submat(2:2:end,1:2:end) = K21; 
    submat(2:2:end,2:2:end) = K22; 

% the kernel for the supported plate integral equation that has to be 
% discretized using smooth quadrature to avoid close evaluations
case {'supported_plate_smooth'}
    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    srcd2 = srcinfo.d2;
    targnorm = targinfo.n;
    targtang = targinfo.d;
    coefs = varargin{1};
    nu = coefs(1);
    
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    
    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
    
    dx = repmat(srctang(1,:),nt,1);
    dy = repmat(srctang(2,:),nt,1);
    
    ds = sqrt(dx.*dx+dy.*dy);
    
    taux = dx./ds;
    tauy = dy./ds;
    
    rx = targ(1,:).' - src(1,:);
    ry = targ(2,:).' - src(2,:);
    r2 = rx.^2 + ry.^2;
    
    dx1 = repmat((targtang(1,:)).',1,ns);
    dy1 = repmat((targtang(2,:)).',1,ns);
    
    ds1 = sqrt(dx1.*dx1+dy1.*dy1); 
    
    tauxtarg = dx1./ds1;
    tauytarg = dy1./ds1;
    
    dx1 = repmat(srctang(1,:),nt,1);
    dy1 = repmat(srctang(2,:),nt,1);
    
    d2x1 = repmat(srcd2(1,:),nt,1);
    d2y1 = repmat(srcd2(2,:),nt,1);
    
    denom = sqrt(dx1.^2+dy1.^2).^3;
    numer = dx1.*d2y1-d2x1.*dy1;
    
    kappa = numer./denom;   
    
    kp = repmat(srcinfo.data(1,:),nt,1);
    kpp = repmat(srcinfo.data(2,:),nt,1);
    
    a1 = 2-nu;
    a2 = (-1+nu)*(7+nu)/(3 - nu);
    a3 = (1-nu)*(3+nu)/(1+nu);
        
    [~, ~, ~, thirdbh, forthbh, fifthbh] = chnk.flex2d.bhgreen(src,targ);
    third = 2*zk^2*thirdbh;
    forth = 2*zk^2*forthbh;
    fifth = 2*zk^2*fifthbh;
    
    submat = -1./(2*zk^2).*(fifth(:,:,1).*nxtarg.^2.*nx.^3 + ...
          fifth(:,:,2).*(2*nxtarg.*nytarg.*nx.^3+3*nxtarg.^2.*nx.^2.*ny) + ...
          fifth(:,:,3).*(nytarg.^2.*nx.^3 + 6*nxtarg.*nytarg.*nx.^2.*ny+3.*nxtarg.^2.*nx.*ny.^2) + ...
          fifth(:,:,4).*(3*nytarg.^2.*nx.^2.*ny+6.*nxtarg.*nytarg.*nx.*ny.^2+nxtarg.^2.*ny.^3)+ ...
          fifth(:,:,5).*(nytarg.*ny.^2.*(3*nytarg.*nx+2*nxtarg.*ny))+ ...
          fifth(:,:,6).*nytarg.^2.*ny.^3) + ... % G_{nx nx ny ny ny}
         -nu./(2*zk^2).*(fifth(:,:,1).*tauxtarg.^2.*nx.^3 + ...
          fifth(:,:,2).*(2*tauxtarg.*tauytarg.*nx.^3+3*tauxtarg.^2.*nx.^2.*ny) + ...
          fifth(:,:,3).*(tauytarg.^2.*nx.^3 + 6*tauxtarg.*tauytarg.*nx.^2.*ny+3.*tauxtarg.^2.*nx.*ny.^2) + ...
          fifth(:,:,4).*(3*tauytarg.^2.*nx.^2.*ny+6.*tauxtarg.*tauytarg.*nx.*ny.^2+tauxtarg.^2.*ny.^3)+ ...
          fifth(:,:,5).*(tauytarg.*ny.^2.*(3*tauytarg.*nx+2*tauxtarg.*ny))+ ...
          fifth(:,:,6).*tauytarg.^2.*ny.^3) + ... % nu*G_{taux taux ny ny ny}
        -a1./(2*zk^2).*(fifth(:,:,1).*nxtarg.^2.*nx.*taux.^2 + ...
          fifth(:,:,2).*(2.*nx.*nxtarg.^2.*taux.*tauy + 2.*nx.*nxtarg.*nytarg.*taux.^2 + nxtarg.^2.*ny.*taux.^2) + ...
          fifth(:,:,3).*(nx.*nxtarg.^2.*tauy.^2 + 4*nx.*nxtarg.*nytarg.*taux.*tauy + 2*nxtarg.^2.*ny.*taux.*tauy + nx.*nytarg.^2.*taux.^2 + 2*nxtarg.*ny.*nytarg.*taux.^2) + ...
          fifth(:,:,4).*(2*nx.*nxtarg.*nytarg.*tauy.^2 + nxtarg.^2.*ny.*tauy.^2 + 2*nx.*nytarg.^2.*taux.*tauy + 4*nxtarg.*ny.*nytarg.*taux.*tauy + ny.*nytarg.^2.*taux.^2) + ...
          fifth(:,:,5).*(nytarg.*tauy.*(2*ny.*nytarg.*taux + 2*nxtarg.*ny.*tauy + nx.*nytarg.*tauy)) + ...
          fifth(:,:,6).*ny.*nytarg.^2.*tauy.^2) + ...  % G_{nx nx ny tauy tauy}
        -nu*a1./(2*zk^2).*(fifth(:,:,1).*tauxtarg.^2.*nx.*taux.^2 + ...
          fifth(:,:,2).*(ny.*taux.^2.*tauxtarg.^2 + 2.*nx.*taux.*tauxtarg.^2.*tauy +2.*nx.*taux.^2.*tauxtarg.*tauytarg) + ...
          fifth(:,:,3).*(nx.*tauxtarg.^2.*tauy.^2 + 4*nx.*tauxtarg.*tauytarg.*taux.*tauy + 2*tauxtarg.^2.*ny.*taux.*tauy + nx.*tauytarg.^2.*taux.^2 + 2*tauxtarg.*ny.*tauytarg.*taux.^2) + ...
          fifth(:,:,4).*(2*nx.*tauxtarg.*tauytarg.*tauy.^2 + tauxtarg.^2.*ny.*tauy.^2 + 2*nx.*tauytarg.^2.*taux.*tauy + 4*tauxtarg.*ny.*tauytarg.*taux.*tauy + ny.*tauytarg.^2.*taux.^2) + ...
          fifth(:,:,5).*(tauytarg.*tauy.*(2*ny.*tauytarg.*taux + 2*tauxtarg.*ny.*tauy + nx.*tauytarg.*tauy)) + ...
          fifth(:,:,6).*ny.*tauytarg.^2.*tauy.^2) + ... % nu*G_{taux taux ny tauy tauy}
         a2.*kappa./(2*zk^2).*(forth(:, :, 1).*(nxtarg.*nxtarg.*nx.*nx) + ...
          forth(:, :, 2).*(2*nx.^2.*nxtarg.*nytarg + 2*nx.*nxtarg.^2.*ny) + ...
          forth(:, :, 3).*(nxtarg.^2.*ny.^2 + 4*nx.*nxtarg.*ny.*nytarg + nx.^2.*nytarg.^2) +...
          forth(:, :, 4).*(2*ny.*nytarg.*(nxtarg.*ny + nx.*nytarg))+...
          forth(:, :, 5).*(nytarg.*nytarg.*ny.*ny)) + ... % kappa*G_{nx nx ny ny}
         nu*a2.*kappa./(2*zk^2).*(forth(:, :, 1).*(tauxtarg.*tauxtarg.*nx.*nx) + ...
          forth(:, :, 2).*(2*nx.^2.*tauxtarg.*tauytarg + 2*nx.*tauxtarg.^2.*ny) + ...
          forth(:, :, 3).*(tauxtarg.^2.*ny.^2 + 4*nx.*tauxtarg.*ny.*tauytarg + nx.^2.*tauytarg.^2) +...
          forth(:, :, 4).*(2*ny.*tauytarg.*(tauxtarg.*ny + nx.*tauytarg))+...
          forth(:, :, 5).*(tauytarg.*tauytarg.*ny.*ny)) + ... % kappa*nu*G_{nx nx ny ny}
        -a3.*kp./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
            third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
            third(:, :, 4).*(nytarg.*nytarg.*tauy)) + ... % kp*G_{nx nx tauy}
         -nu*a3.*kp./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
            third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
            third(:, :, 4).*(tauytarg.*tauytarg.*tauy)) ;  % kp*nu*G_{taux taux tauy}

    submat(r2 == 0) = (nu - 1)*(12*kappa(r2 == 0).^3*(nu^2 - nu + 4) + kpp(r2 == 0)*(-5*nu^2 + 4*nu + 33))/(48*pi*(nu - 3)) ; % diagonal replacement

% supported plate kernels for plotting
case {'supported_plate_eval'}
    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    srcd2 = srcinfo.d2;
    coefs = varargin{1};
    nu = coefs(1);
    
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    
    dx = repmat(srctang(1,:),nt,1);
    dy = repmat(srctang(2,:),nt,1);
    
    ds = sqrt(dx.*dx+dy.*dy);
    
    taux = dx./ds;
    tauy = dy./ds;
    
    dx1 = repmat(srctang(1,:),nt,1);
    dy1 = repmat(srctang(2,:),nt,1);
    
    d2x1 = repmat(srcd2(1,:),nt,1);
    d2y1 = repmat(srcd2(2,:),nt,1);
    
    denom = sqrt(dx1.^2+dy1.^2).^3;
    numer = dx1.*d2y1-d2x1.*dy1;
    
    kappa = numer./denom; 
    
    kp = repmat(srcinfo.data(1,:),nt,1);
    
    a1 = 2-nu;
    a2 = (-1+nu)*(7+nu)/(3 - nu);
    a3 = (1-nu)*(3+nu)/(1+nu);
    
    [~, grad, hess, third] = chnk.flex2d.hkdiffgreen(zk, src, targ, false); 
    
    K1 = -1/(2*zk^2)*(third(:,:,1).*nx.^3 + 3*third(:,:,2).*nx.^2.*ny + 3*third(:,:,3).*nx.*ny.^2 + third(:,:,4).*ny.^3) + ...
         -a1/(2*zk^2)*(third(:,:,1).*nx.*taux.^2 + third(:,:,2).*(ny.*taux.^2 + 2*nx.*taux.*tauy) + third(:,:,3).*(nx.*tauy.^2 + 2*ny.*taux.*tauy) + third(:,:,4).*ny.*tauy.^2) + ...
         a2*kappa./(2*zk^2).*(hess(:,:,1).*nx.^2 + 2*hess(:,:,2).*nx.*ny + hess(:,:,3).*ny.^2) + ...
         -a3*kp./(2*zk^2).*(grad(:,:,1).*taux + grad(:,:,2).*tauy);

    K2 = -1/(2*zk^2).*(grad(:,:,1).*nx + grad(:,:,2).*ny);

    submat = zeros(nt,2*ns);

    submat(:,1:2:end) = K1;
    submat(:,2:2:end) = K2;

end

end


