function submat= kern(zk,srcinfo,targinfo,type,kappa,d,Sn,s0_l,sn_l,l,ising,varargin)

% 
% see also CHNK.FLEX2D.KERN
  
src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);
nkappa = length(kappa);

%%% STANDARD LAYER POTENTIALS

switch lower(type)
case {'s', 'single'} % flexural wave single layer

   val = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,0);  
   submat = 1/(2*zk^2).*val;

case {'sp', 'sprime'} % normal derivative of flexural wave single layer

    targnorm = targinfo.n;
    nxtarg = repmat(reshape(targnorm(1,:),1,nt,1),nkappa,1,ns);
    nxtarg = reshape(nxtarg,[],ns);
    nytarg = repmat(reshape(targnorm(2,:),1,nt,1),nkappa,1,ns);
    nytarg = reshape(nytarg,[],ns);

    [~,grad] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l);  
    submat = 1/(2*zk^2).*(grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);

case {'d', 'double'} % normal derivative of flexural wave single layer

    srcnorm = srcinfo.n;
    nx = repmat(srcnorm(1,:),nkappa*nt,1);
    ny = repmat(srcnorm(2,:),nkappa*nt,1);
    
    [~,grad] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l);  
    submat = -1/(2*zk^2).*(grad(:,:,1).*nx + grad(:,:,2).*ny);


%%% CLAMPED PLATE KERNELS

% boundary conditions applied to a point source
case {'clamped_plate_bcs'}
    targnorm = targinfo.n;
    nxtarg = repmat(reshape(targnorm(1,:),1,nt,1),nkappa,1,ns);
    nxtarg = reshape(nxtarg,[],ns);
    nytarg = repmat(reshape(targnorm(2,:),1,nt,1),nkappa,1,ns);
    nytarg = reshape(nytarg,[],ns);
    
    [val, grad] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,0);  
    
    firstbc = 1/(2*zk^2).*val ;
    secondbc = 1/(2*zk^2).*(grad(:, :, 1).*nxtarg + grad(:, :, 2).*nytarg);
   
    submat = zeros(nkappa,2,nt,ns);
    submat(:,1,:,:) = reshape(firstbc,nkappa,1,nt,[]);
    submat(:,2,:,:) = reshape(secondbc,nkappa,1,nt,[]);
    submat = reshape(submat, [],ns);

case {'clamped_plate_bcs_trx'}
    submat = zeros(nkappa,2,nt,4*ns);

    targnorm = targinfo.n;
    nxtarg = repmat(reshape(targnorm(1,:),1,nt,1),nkappa,1,ns);
    nxtarg = reshape(nxtarg,[],ns);
    nytarg = repmat(reshape(targnorm(2,:),1,nt,1),nkappa,1,ns);
    nytarg = reshape(nytarg,[],ns);
    
    [val, grad, hess, third,fourth] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,1);  
    
    firstbc = 1/(2*zk^2).*val ;
    secondbc = 1/(2*zk^2).*(grad(:, :, 1).*nxtarg + grad(:, :, 2).*nytarg);
   
    submat(:,1,:,4:4:end) = reshape(firstbc,nkappa,1,nt,[]);
    submat(:,2,:,4:4:end) = reshape(secondbc,nkappa,1,nt,[]);

    firstbc = 1/(2*zk^2).*grad(:, :, 1) ;
    secondbc = 1/(2*zk^2).*(hess(:, :, 1).*nxtarg + hess(:, :, 2).*nytarg);
   
    submat(:,1,:,3:4:end) = -reshape(firstbc,nkappa,1,nt,[]);
    submat(:,2,:,3:4:end) = -reshape(secondbc,nkappa,1,nt,[]);

    firstbc = 1/(2*zk^2).*hess(:, :, 1) ;
    secondbc = 1/(2*zk^2).*(third(:, :, 1).*nxtarg + third(:, :, 2).*nytarg);

    firstbc = firstbc + 2*1/(2*zk^2).*hess(:, :, 3) ;
    secondbc = secondbc + 2*(1/(2*zk^2).*(third(:, :, 3).*nxtarg + third(:, :, 4).*nytarg));

    submat(:,1,:,2:4:end) = reshape(firstbc,nkappa,1,nt,[]);
    submat(:,2,:,2:4:end) = reshape(secondbc,nkappa,1,nt,[]);

    firstbc = 1/(2*zk^2).*third(:, :, 1) ;
    secondbc = 1/(2*zk^2).*(fourth(:, :, 1).*nxtarg + fourth(:, :, 2).*nytarg);

    firstbc = firstbc+2*1/(2*zk^2).*third(:, :, 3) ;
    secondbc = secondbc + 2*(1/(2*zk^2).*(fourth(:, :, 3).*nxtarg + fourth(:, :, 4).*nytarg));


    submat(:,1,:,1:4:end) = -reshape(firstbc,nkappa,1,nt,[]);
    submat(:,2,:,1:4:end) = -reshape(secondbc,nkappa,1,nt,[]);
    
    submat = reshape(submat, [],4*ns);

% kernels for the clamped plate integral equation
case {'clamped_plate'}
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;

    nx = repmat(srcnorm(1,:),nkappa*nt,1);
    ny = repmat(srcnorm(2,:),nkappa*nt,1);
   
    targnorm = targinfo.n;
    nxtarg = repmat(reshape(targnorm(1,:),1,nt,1),nkappa,1,ns);
    nxtarg = reshape(nxtarg,[],ns);
    nytarg = repmat(reshape(targnorm(2,:),1,nt,1),nkappa,1,ns);
    nytarg = reshape(nytarg,[],ns);

    if length(varargin) > 0
        nsub = varargin{1};
    else 
        nsub = 0;
    end

   
   [~, ~, hess, third, fourth] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,0,nsub);  
   
   dx = repmat(srctang(1,:),nkappa*nt,1);
   dy = repmat(srctang(2,:),nkappa*nt,1);
    
   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds;
   tauy = dy./ds;
   
   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   % r2 = rx.^2 + ry.^2;
   % 
   % rn = rx.*nx + ry.*ny;
   % rtau = rx.*taux + ry.*tauy;
   % ntargtau = nxtarg.*taux + nytarg.*tauy;
   % 
   % rntarg = rx.*nxtarg + ry.*nytarg;

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
          fourth(:, :, 5).*(ny.*tauy.*tauy.*nytarg)));

   K22 = -(1/(2*zk^2).*(third(:,:, 1).*(nx.*nx.*nxtarg) +third(:, :, 2).*(nx.*nx.*nytarg + 2*nx.*ny.*nxtarg) + third(:, :, 3).*(2*nx.*ny.*nytarg + ny.*ny.*nxtarg)+...
         third(:, :,4).*(ny.*ny.*nytarg))) + ...
         (1/(2*zk^2).*(third(:,:, 1).*(taux.*taux.*nxtarg) +third(:, :, 2).*(taux.*taux.*nytarg + 2*taux.*tauy.*nxtarg) + third(:, :, 3).*(2*taux.*tauy.*nytarg + tauy.*tauy.*nxtarg)+...
         third(:, :,4).*(tauy.*tauy.*nytarg)));

  submat = zeros(nkappa,2,nt,2*ns);
  submat(:,1,:,1:2:2*ns) = reshape(K11,nkappa,1,nt,[]);
  submat(:,1,:,2:2:2*ns) = reshape(K12,nkappa,1,nt,[]);
  submat(:,2,:,1:2:2*ns) = reshape(K21,nkappa,1,nt,[]);
  submat(:,2,:,2:2:2*ns) = reshape(K22,nkappa,1,nt,[]);
  submat = reshape(submat, [],2*ns);

  % submat = zeros(2*nt,2*ns);
  % 
  % submat(1:2:end,1:2:end) = K11;
  % submat(1:2:end,2:2:end) = K12;
  % 
  % submat(2:2:end,1:2:end) = K21;
  % submat(2:2:end,2:2:end) = K22;

% clamped plate kernels for plotting
case {'clamped_plate_eval'}

    submat = zeros(nkappa*nt,2*ns);

    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    nx = repmat(srcnorm(1,:),nkappa*nt,1);
    ny = repmat(srcnorm(2,:),nkappa*nt,1);
   dx = repmat(srctang(1,:),nkappa*nt,1);
   dy = repmat(srctang(2,:),nkappa*nt,1);
    ds = sqrt(dx.*dx+dy.*dy);

    taux = dx./ds;
    tauy = dy./ds;

    [~, ~, hess, third] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,0);  

    K1 = -(1/(2*zk^2).*(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
       third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny)) ) - ...
       (3/(2*zk^2).*(third(:, :, 1).*(nx.*taux.*taux) + third(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       third(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + third(:, :, 4).*(ny.*tauy.*tauy)));  % G_{ny ny ny} + 3G_{ny tauy tauy}

    K2 =  -(1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny)))+...
          (1/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy))); % -G_{ny ny}  + G_{tauy tauy}

    submat(:,1:2:end) = K1;
    submat(:,2:2:end) = K2;


case {'clamped_plate_eval_trx'}

    submat = zeros(nkappa,4,nt,2*ns);

    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    nx = repmat(srcnorm(1,:),nkappa*nt,1);
    ny = repmat(srcnorm(2,:),nkappa*nt,1);
   dx = repmat(srctang(1,:),nkappa*nt,1);
   dy = repmat(srctang(2,:),nkappa*nt,1);
    ds = sqrt(dx.*dx+dy.*dy);

    taux = dx./ds;
    tauy = dy./ds;

    [~, ~, hess, third,fourth,fifth,sixth] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,1);  

    K1 = -(1/(2*zk^2).*(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
       third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny)) ) - ...
       (3/(2*zk^2).*(third(:, :, 1).*(nx.*taux.*taux) + third(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       third(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + third(:, :, 4).*(ny.*tauy.*tauy)));  % G_{ny ny ny} + 3G_{ny tauy tauy}

    K2 =  -(1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny)))+...
          (1/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy))); % -G_{ny ny}  + G_{tauy tauy}

    submat(:,1,:,1:2:end) = reshape(K1,nkappa,1,nt,[]);
    submat(:,1,:,2:2:end) = reshape(K2,nkappa,1,nt,[]);

    K1 = -(1/(2*zk^2).*(fourth(:, :, 1).*(nx.*nx.*nx) + fourth(:, :, 2).*(3*nx.*nx.*ny) +...
       fourth(:, :, 3).*(3*nx.*ny.*ny) + fourth(:, :, 4).*(ny.*ny.*ny)) ) - ...
       (3/(2*zk^2).*(fourth(:, :, 1).*(nx.*taux.*taux) + fourth(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       fourth(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + fourth(:, :, 4).*(ny.*tauy.*tauy)));  % G_{ny ny ny} + 3G_{ny tauy tauy}

    K2 =  -(1/(2*zk^2).*(third(:, :, 1).*(nx.*nx) + third(:, :, 2).*(2*nx.*ny) + third(:, :, 3).*(ny.*ny)))+...
          (1/(2*zk^2).*(third(:, :, 1).*(taux.*taux) + third(:, :, 2).*(2*taux.*tauy) + third(:, :, 3).*(tauy.*tauy))); % -G_{ny ny}  + G_{tauy tauy}

    submat(:,2,:,1:2:end) = reshape(K1,nkappa,1,nt,[]);
    submat(:,2,:,2:2:end) = reshape(K2,nkappa,1,nt,[]);

    K1 = -(1/(2*zk^2).*(fifth(:, :, 1).*(nx.*nx.*nx) + fifth(:, :, 2).*(3*nx.*nx.*ny) +...
       fifth(:, :, 3).*(3*nx.*ny.*ny) + fifth(:, :, 4).*(ny.*ny.*ny)) ) - ...
       (3/(2*zk^2).*(fifth(:, :, 1).*(nx.*taux.*taux) + fifth(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       fifth(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + fifth(:, :, 4).*(ny.*tauy.*tauy)));  % G_{ny ny ny} + 3G_{ny tauy tauy}

    K2 =  -(1/(2*zk^2).*(fourth(:, :, 1).*(nx.*nx) + fourth(:, :, 2).*(2*nx.*ny) + fourth(:, :, 3).*(ny.*ny)))+...
          (1/(2*zk^2).*(fourth(:, :, 1).*(taux.*taux) + fourth(:, :, 2).*(2*taux.*tauy) + fourth(:, :, 3).*(tauy.*tauy))); % -G_{ny ny}  + G_{tauy tauy}

    submat(:,3,:,1:2:end) = reshape(K1,nkappa,1,nt,[]);
    submat(:,3,:,2:2:end) = reshape(K2,nkappa,1,nt,[]);

    K1 = -(1/(2*zk^2).*(sixth(:, :, 1).*(nx.*nx.*nx) + sixth(:, :, 2).*(3*nx.*nx.*ny) +...
       sixth(:, :, 3).*(3*nx.*ny.*ny) + sixth(:, :, 4).*(ny.*ny.*ny)) ) - ...
       (3/(2*zk^2).*(sixth(:, :, 1).*(nx.*taux.*taux) + sixth(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       sixth(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + sixth(:, :, 4).*(ny.*tauy.*tauy)));  % G_{ny ny ny} + 3G_{ny tauy tauy}

    K2 =  -(1/(2*zk^2).*(fifth(:, :, 1).*(nx.*nx) + fifth(:, :, 2).*(2*nx.*ny) + fifth(:, :, 3).*(ny.*ny)))+...
          (1/(2*zk^2).*(fifth(:, :, 1).*(taux.*taux) + fifth(:, :, 2).*(2*taux.*tauy) + fifth(:, :, 3).*(tauy.*tauy))); % -G_{ny ny}  + G_{tauy tauy}

    submat(:,4,:,1:2:end) = reshape(K1,nkappa,1,nt,[]);
    submat(:,4,:,2:2:end) = reshape(K2,nkappa,1,nt,[]);

    submat = reshape(submat, [],2*ns);


%%% FREE PLATE KERNELS

% boundary conditions applied to a point source
case {'free_plate_bcs'}
    targtang = targinfo.d;
    targd2 = targinfo.d2;
    nu = varargin{1};
    
    [~, ~, hess, third] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,0);  

    targnorm = targinfo.n;
    nxtarg = repmat(reshape(targnorm(1,:),1,nt,1),nkappa,1,ns);
    nxtarg = reshape(nxtarg,[],ns);
    nytarg = repmat(reshape(targnorm(2,:),1,nt,1),nkappa,1,ns);
    nytarg = reshape(nytarg,[],ns);

    dx1 = repmat(reshape(targtang(1,:),1,nt,1),nkappa,1,ns);
    dx1 = reshape(dx1,[],ns);
    dy1 = repmat(reshape(targtang(2,:),1,nt,1),nkappa,1,ns);
    dy1 = reshape(dy1,[],ns);

    % dx1 = repmat((targtang(1,:)).',1,ns);
    % dy1 = repmat((targtang(2,:)).',1,ns);
    
    ds1 = sqrt(dx1.*dx1+dy1.*dy1); 
    
    % d2x1 = repmat((targd2(1,:)).',1,ns);
    % d2y1 = repmat((targd2(2,:)).',1,ns);
    
    d2x1 = repmat(reshape(targd2(1,:),1,nt,1),nkappa,1,ns);
    d2x1 = reshape(d2x1,[],ns);
    d2y1 = repmat(reshape(targd2(2,:),1,nt,1),nkappa,1,ns);
    d2y1 = reshape(d2y1,[],ns);

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

    submat = zeros(nkappa,2,nt,ns);
    submat(:,1,:,:) = reshape(firstbc,nkappa,1,nt,[]);
    submat(:,2,:,:) = reshape(secondbc,nkappa,1,nt,[]);
    submat = reshape(submat, [],ns);

case {'free_plate_bcs_trx'}
    % transmission boundary to object
    submat = zeros(nkappa,2,nt,4*ns);

    targtang = targinfo.d;
    targd2 = targinfo.d2;
    nu = varargin{1};
    
    [~, ~, hess, third, fourth, fifth, sixth] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,1);  

    targnorm = targinfo.n;
    nxtarg = repmat(reshape(targnorm(1,:),1,nt,1),nkappa,1,ns);
    nxtarg = reshape(nxtarg,[],ns);
    nytarg = repmat(reshape(targnorm(2,:),1,nt,1),nkappa,1,ns);
    nytarg = reshape(nytarg,[],ns);

    dx1 = repmat(reshape(targtang(1,:),1,nt,1),nkappa,1,ns);
    dx1 = reshape(dx1,[],ns);
    dy1 = repmat(reshape(targtang(2,:),1,nt,1),nkappa,1,ns);
    dy1 = reshape(dy1,[],ns);

    % dx1 = repmat((targtang(1,:)).',1,ns);
    % dy1 = repmat((targtang(2,:)).',1,ns);
    
    ds1 = sqrt(dx1.*dx1+dy1.*dy1); 
    
    % d2x1 = repmat((targd2(1,:)).',1,ns);
    % d2y1 = repmat((targd2(2,:)).',1,ns);
    
    d2x1 = repmat(reshape(targd2(1,:),1,nt,1),nkappa,1,ns);
    d2x1 = reshape(d2x1,[],ns);
    d2y1 = repmat(reshape(targd2(2,:),1,nt,1),nkappa,1,ns);
    d2y1 = reshape(d2y1,[],ns);

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

    submat(:,1,:,4:4:end) = reshape(firstbc,nkappa,1,nt,[]);
    submat(:,2,:,4:4:end) = reshape(secondbc,nkappa,1,nt,[]);

    firstbc = 1/(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg) + third(:, :, 2).*(2*nxtarg.*nytarg) + third(:, :, 3).*(nytarg.*nytarg))+...
    nu/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg) + third(:, :, 2).*(2*tauxtarg.*tauytarg) + third(:, :, 3).*(tauytarg.*tauytarg));
    
    secondbc = 1./(2*zk^2).*(fourth(:, :, 1).*(nxtarg.*nxtarg.*nxtarg) + fourth(:, :, 2).*(3*nxtarg.*nxtarg.*nytarg) +...
    fourth(:, :, 3).*(3*nxtarg.*nytarg.*nytarg) + fourth(:, :, 4).*(nytarg.*nytarg.*nytarg))+...
    (2-nu)/(2*zk^2).*(fourth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg) + fourth(:, :, 2).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
    fourth(:, :, 3).*(2*tauxtarg.*tauytarg.*nytarg+ tauytarg.*tauytarg.*nxtarg) +...
    + fourth(:, :, 4).*(tauytarg.*tauytarg.*nytarg))+...
    (1-nu).*kappatarg.*(1/(2*zk^2).*(third(:, :, 1).*tauxtarg.*tauxtarg + third(:, :, 2).*(2*tauxtarg.*tauytarg) + third(:, :, 3).*tauytarg.*tauytarg)-...
    (1/(2*zk^2).*(third(:, :, 1).*nxtarg.*nxtarg + third(:, :, 2).*(2*nxtarg.*nytarg) + third(:, :, 3).*nytarg.*nytarg)));

    submat(:,1,:,3:4:end) = -reshape(firstbc,nkappa,1,nt,[]);
    submat(:,2,:,3:4:end) = -reshape(secondbc,nkappa,1,nt,[]);

    firstbc = 1/(2*zk^2).*(fourth(:, :, 1).*(nxtarg.*nxtarg) + fourth(:, :, 2).*(2*nxtarg.*nytarg) + fourth(:, :, 3).*(nytarg.*nytarg))+...
    nu/(2*zk^2).*(fourth(:, :, 1).*(tauxtarg.*tauxtarg) + fourth(:, :, 2).*(2*tauxtarg.*tauytarg) + fourth(:, :, 3).*(tauytarg.*tauytarg));
    
    firstbc = firstbc + 2*(1/(2*zk^2).*(fourth(:, :, 3).*(nxtarg.*nxtarg) + fourth(:, :, 4).*(2*nxtarg.*nytarg) + fourth(:, :, 5).*(nytarg.*nytarg))+...
    nu/(2*zk^2).*(fourth(:, :, 3).*(tauxtarg.*tauxtarg) + fourth(:, :, 4).*(2*tauxtarg.*tauytarg) + fourth(:, :, 5).*(tauytarg.*tauytarg)));
    
    secondbc = 1./(2*zk^2).*(fifth(:, :, 1).*(nxtarg.*nxtarg.*nxtarg) + fifth(:, :, 2).*(3*nxtarg.*nxtarg.*nytarg) +...
    fifth(:, :, 3).*(3*nxtarg.*nytarg.*nytarg) + fifth(:, :, 4).*(nytarg.*nytarg.*nytarg))+...
    (2-nu)/(2*zk^2).*(fifth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg) + fifth(:, :, 2).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
    fifth(:, :, 3).*(2*tauxtarg.*tauytarg.*nytarg+ tauytarg.*tauytarg.*nxtarg) +...
    + fifth(:, :, 4).*(tauytarg.*tauytarg.*nytarg))+...
    (1-nu).*kappatarg.*(1/(2*zk^2).*(fourth(:, :, 1).*tauxtarg.*tauxtarg + fourth(:, :, 2).*(2*tauxtarg.*tauytarg) + fourth(:, :, 3).*tauytarg.*tauytarg)-...
    (1/(2*zk^2).*(fourth(:, :, 1).*nxtarg.*nxtarg + fourth(:, :, 2).*(2*nxtarg.*nytarg) + fourth(:, :, 3).*nytarg.*nytarg)));

    secondbc = secondbc + 2*(1./(2*zk^2).*(fifth(:, :, 3).*(nxtarg.*nxtarg.*nxtarg) + fifth(:, :, 4).*(3*nxtarg.*nxtarg.*nytarg) +...
    fifth(:, :, 5).*(3*nxtarg.*nytarg.*nytarg) + fifth(:, :, 6).*(nytarg.*nytarg.*nytarg))+...
    (2-nu)/(2*zk^2).*(fifth(:, :, 3).*(tauxtarg.*tauxtarg.*nxtarg) + fifth(:, :, 4).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
    fifth(:, :, 5).*(2*tauxtarg.*tauytarg.*nytarg+ tauytarg.*tauytarg.*nxtarg) +...
    + fifth(:, :, 6).*(tauytarg.*tauytarg.*nytarg))+...
    (1-nu).*kappatarg.*(1/(2*zk^2).*(fourth(:, :, 3).*tauxtarg.*tauxtarg + fourth(:, :, 4).*(2*tauxtarg.*tauytarg) + fourth(:, :, 5).*tauytarg.*tauytarg)-...
    (1/(2*zk^2).*(fourth(:, :, 3).*nxtarg.*nxtarg + fourth(:, :, 4).*(2*nxtarg.*nytarg) + fourth(:, :, 5).*nytarg.*nytarg))));


    submat(:,1,:,2:4:end) = reshape(firstbc,nkappa,1,nt,[]);
    submat(:,2,:,2:4:end) = reshape(secondbc,nkappa,1,nt,[]);

    firstbc = 1/(2*zk^2).*(fifth(:, :, 1).*(nxtarg.*nxtarg) + fifth(:, :, 2).*(2*nxtarg.*nytarg) + fifth(:, :, 3).*(nytarg.*nytarg))+...
    nu/(2*zk^2).*(fifth(:, :, 1).*(tauxtarg.*tauxtarg) + fifth(:, :, 2).*(2*tauxtarg.*tauytarg) + fifth(:, :, 3).*(tauytarg.*tauytarg));
    
    firstbc = firstbc + 2*(1/(2*zk^2).*(fifth(:, :, 3).*(nxtarg.*nxtarg) + fifth(:, :, 4).*(2*nxtarg.*nytarg) + fifth(:, :, 5).*(nytarg.*nytarg))+...
    nu/(2*zk^2).*(fifth(:, :, 3).*(tauxtarg.*tauxtarg) + fifth(:, :, 4).*(2*tauxtarg.*tauytarg) + fifth(:, :, 5).*(tauytarg.*tauytarg)));
    
    secondbc = 1./(2*zk^2).*(sixth(:, :, 1).*(nxtarg.*nxtarg.*nxtarg) + sixth(:, :, 2).*(3*nxtarg.*nxtarg.*nytarg) +...
    sixth(:, :, 3).*(3*nxtarg.*nytarg.*nytarg) + sixth(:, :, 4).*(nytarg.*nytarg.*nytarg))+...
    (2-nu)/(2*zk^2).*(sixth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg) + sixth(:, :, 2).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
    sixth(:, :, 3).*(2*tauxtarg.*tauytarg.*nytarg+ tauytarg.*tauytarg.*nxtarg) +...
    + sixth(:, :, 4).*(tauytarg.*tauytarg.*nytarg))+...
    (1-nu).*kappatarg.*(1/(2*zk^2).*(fifth(:, :, 1).*tauxtarg.*tauxtarg + fifth(:, :, 2).*(2*tauxtarg.*tauytarg) + fifth(:, :, 3).*tauytarg.*tauytarg)-...
    (1/(2*zk^2).*(fifth(:, :, 1).*nxtarg.*nxtarg + fifth(:, :, 2).*(2*nxtarg.*nytarg) + fifth(:, :, 3).*nytarg.*nytarg)));

    secondbc = secondbc+2*(1./(2*zk^2).*(sixth(:, :, 3).*(nxtarg.*nxtarg.*nxtarg) + sixth(:, :, 4).*(3*nxtarg.*nxtarg.*nytarg) +...
    sixth(:, :, 5).*(3*nxtarg.*nytarg.*nytarg) + sixth(:, :, 6).*(nytarg.*nytarg.*nytarg))+...
    (2-nu)/(2*zk^2).*(sixth(:, :, 3).*(tauxtarg.*tauxtarg.*nxtarg) + sixth(:, :, 4).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
    sixth(:, :, 5).*(2*tauxtarg.*tauytarg.*nytarg+ tauytarg.*tauytarg.*nxtarg) +...
    + sixth(:, :, 6).*(tauytarg.*tauytarg.*nytarg))+...
    (1-nu).*kappatarg.*(1/(2*zk^2).*(fifth(:, :, 3).*tauxtarg.*tauxtarg + fifth(:, :, 4).*(2*tauxtarg.*tauytarg) + fifth(:, :, 5).*tauytarg.*tauytarg)-...
    (1/(2*zk^2).*(fifth(:, :, 3).*nxtarg.*nxtarg + fifth(:, :, 4).*(2*nxtarg.*nytarg) + fifth(:, :, 5).*nytarg.*nytarg))));


    submat(:,1,:,1:4:end) = -reshape(firstbc,nkappa,1,nt,[]);
    submat(:,2,:,1:4:end) = -reshape(secondbc,nkappa,1,nt,[]);

    submat = reshape(submat, [],4*ns);

% kernels for the free plate integral equation 
case {'free_plate'}
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   targd2 = targinfo.d2;
   nu = varargin{1};

   if length(varargin) > 1
        nsub = varargin{2};
   else
        nsub = 0;
   end

   [~, ~, hess, third, fourth] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,0,nsub);  

    nx = repmat(srcnorm(1,:),nkappa*nt,1);
    ny = repmat(srcnorm(2,:),nkappa*nt,1);

    nxtarg = repmat(reshape(targnorm(1,:),1,nt,1),nkappa,1,ns);
    nxtarg = reshape(nxtarg,[],ns);
    nytarg = repmat(reshape(targnorm(2,:),1,nt,1),nkappa,1,ns);
    nytarg = reshape(nytarg,[],ns);

   dx = repmat(srctang(1,:),nkappa*nt,1);
   dy = repmat(srctang(2,:),nkappa*nt,1);

   ds = sqrt(dx.*dx+dy.*dy); 

   taux = dx ./ ds;
   tauy = dy ./ ds;

    dx1 = repmat(reshape(targtang(1,:),1,nt,1),nkappa,1,ns);
    dx1 = reshape(dx1,[],ns);
    dy1 = repmat(reshape(targtang(2,:),1,nt,1),nkappa,1,ns);
    dy1 = reshape(dy1,[],ns);

    % dx1 = repmat((targtang(1,:)).',1,ns);
    % dy1 = repmat((targtang(2,:)).',1,ns);
    
    ds1 = sqrt(dx1.*dx1+dy1.*dy1); 
    
    % d2x1 = repmat((targd2(1,:)).',1,ns);
    % d2y1 = repmat((targd2(2,:)).',1,ns);
    
    d2x1 = repmat(reshape(targd2(1,:),1,nt,1),nkappa,1,ns);
    d2x1 = reshape(d2x1,[],ns);
    d2y1 = repmat(reshape(targd2(2,:),1,nt,1),nkappa,1,ns);
    d2y1 = reshape(d2y1,[],ns);

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   denom = sqrt(dx1.^2 + dy1.^2).^3;
   numer = dx1.*d2y1 - d2x1.*dy1;

   kappatarg = numer ./ denom; % target curvature

   hilb = chnk.lap2dquas.kern(srcinfo,targinfo,'hilb',kappa,d,s0_l,sn_l,l,0,nsub);
   hilbp = chnk.lap2dquas.kern(srcinfo,targinfo,'hilbprime',kappa,d,s0_l,sn_l,l,0,nsub);
   hilb = (1+nu) * hilb / 2;
   hilbp = (1+nu) * hilbp / 2;
   
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
          fourth(:, :, 5).*(tauytarg.*tauytarg.*nytarg.*ny)) ) ...          
          -hilbp/2;

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
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy))) + (1+ nu)/4*hilb;

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
         fourth(:, :, 5).*(tauytarg.*tauytarg.*nytarg.*tauy))) ;
    
    % submat = zeros(4*nkappa*nt,2*ns);
    % 
    % submat(1:4:end,1:2:end) = K11;
    % submat(1:4:end,2:2:end) = K12;
    % 
    % submat(2:4:end,1:2:end) = K21;
    % submat(2:4:end,2:2:end) = K22;
    % 
    % submat(3:4:end,1:2:end) = K11H;
    % submat(4:4:end,1:2:end) = K21H;
      
    submat = zeros(nkappa,4,nt,2*ns);
  submat(:,1,:,1:2:2*ns) = reshape(K11,nkappa,1,nt,[]);
  submat(:,1,:,2:2:2*ns) = reshape(K12,nkappa,1,nt,[]);
  submat(:,2,:,1:2:2*ns) = reshape(K21,nkappa,1,nt,[]);
  submat(:,2,:,2:2:2*ns) = reshape(K22,nkappa,1,nt,[]);
  submat(:,3,:,1:2:2*ns) = reshape(K11H,nkappa,1,nt,[]);
  submat(:,4,:,1:2:2*ns) = reshape(K21H,nkappa,1,nt,[]);
  submat = reshape(submat, [],2*ns);

% free plate kernels used for plotting 
case {'free_plate_eval'}
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   nu = varargin{1};

   [val,grad] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,0);  
    nx = repmat(srcnorm(1,:),nkappa*nt,1);
    ny = repmat(srcnorm(2,:),nkappa*nt,1);

   dx = repmat(srctang(1,:),nkappa*nt,1);
   dy = repmat(srctang(2,:),nkappa*nt,1);

   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds; 
   tauy = dy./ds;
  
   K1 = (-1/(2*zk^2).*(grad(:, :, 1).*(nx) + grad(:, :, 2).*ny)); 
   K1H = ((1 + nu)/2).*(-1/(2*zk^2).*(grad(:, :, 1).*(taux) + grad(:, :, 2).*tauy));                    % G_{tauy}
   K2 = 1/(2*zk^2).*val;

   submat = zeros(nkappa*nt,3*ns);
   submat(:,1:3:end) = K1;
   submat(:,2:3:end) = K1H;
   submat(:,3:3:end) = K2;

case {'free_plate_eval_trx'}

    submat = zeros(nkappa,4,nt,3*ns);

    srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   nu = varargin{1};

   [val,grad,hess,third,fourth] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,1);  
    nx = repmat(srcnorm(1,:),nkappa*nt,1);
    ny = repmat(srcnorm(2,:),nkappa*nt,1);

   dx = repmat(srctang(1,:),nkappa*nt,1);
   dy = repmat(srctang(2,:),nkappa*nt,1);

   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds; 
   tauy = dy./ds;
  
   K1 = (-1/(2*zk^2).*(grad(:, :, 1).*(nx) + grad(:, :, 2).*ny)); 
   K1H = ((1 + nu)/2).*(-1/(2*zk^2).*(grad(:, :, 1).*(taux) + grad(:, :, 2).*tauy));                    % G_{tauy}
   K2 = 1/(2*zk^2).*val;

   submat(:,1,:,1:3:end) = reshape(K1,nkappa,1,nt,[]);
   submat(:,1,:,2:3:end) = reshape(K1H,nkappa,1,nt,[]);
   submat(:,1,:,3:3:end) = reshape(K2,nkappa,1,nt,[]);

   K1 = (-1/(2*zk^2).*(hess(:, :, 1).*(nx) + hess(:, :, 2).*ny)); 
   K1H = ((1 + nu)/2).*(-1/(2*zk^2).*(hess(:, :, 1).*(taux) + hess(:, :, 2).*tauy));                    % G_{tauy}
   K2 = 1/(2*zk^2).*grad(:,:,1);

   submat(:,2,:,1:3:end) = reshape(K1,nkappa,1,nt,[]);
   submat(:,2,:,2:3:end) = reshape(K1H,nkappa,1,nt,[]);
   submat(:,2,:,3:3:end) = reshape(K2,nkappa,1,nt,[]);

   K1 = (-1/(2*zk^2).*(third(:, :, 1).*(nx) + third(:, :, 2).*ny)); 
   K1H = ((1 + nu)/2).*(-1/(2*zk^2).*(third(:, :, 1).*(taux) + third(:, :, 2).*tauy));                    % G_{tauy}
   K2 = 1/(2*zk^2).*hess(:,:,1);

   submat(:,3,:,1:3:end) = reshape(K1,nkappa,1,nt,[]);
   submat(:,3,:,2:3:end) = reshape(K1H,nkappa,1,nt,[]);
   submat(:,3,:,3:3:end) = reshape(K2,nkappa,1,nt,[]);

   K1 = (-1/(2*zk^2).*(fourth(:, :, 1).*(nx) + fourth(:, :, 2).*ny)); 
   K1H = ((1 + nu)/2).*(-1/(2*zk^2).*(fourth(:, :, 1).*(taux) + fourth(:, :, 2).*tauy));                    % G_{tauy}
   K2 = 1/(2*zk^2).*third(:,:,1);

   submat(:,4,:,1:3:end) = reshape(K1,nkappa,1,nt,[]);
   submat(:,4,:,2:3:end) = reshape(K1H,nkappa,1,nt,[]);
   submat(:,4,:,3:3:end) = reshape(K2,nkappa,1,nt,[]);

   submat = reshape(submat, [],3*ns);

% boundary conditions applied to a point source
case {'supported_plate_bcs'}
   targnorm = targinfo.n;
   targtang = targinfo.d;

    nxtarg = repmat(reshape(targnorm(1,:),1,nt,1),nkappa,1,ns);
    nxtarg = reshape(nxtarg,[],ns);
    nytarg = repmat(reshape(targnorm(2,:),1,nt,1),nkappa,1,ns);
    nytarg = reshape(nytarg,[],ns);
    dx = repmat(reshape(targtang(1,:),1,nt,1),nkappa,1,ns);
    dx = reshape(dx,[],ns);
    dy = repmat(reshape(targtang(2,:),1,nt,1),nkappa,1,ns);
    dy = reshape(dy,[],ns);
    % dx = targinfo.d(1,:).';
    % dy = targinfo.d(2,:).';
    ds = sqrt(dx.*dx+dy.*dy);
    tauxtarg = (dx./ds);                                                                       % normalization
    tauytarg = (dy./ds);

    nu = varargin{1};

    [val, ~, hess] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,0);  
    
    firstbc = 1/(2*zk^2).*val ;
    
    secondbc = 1/(2*zk^2).*(hess(:, :, 1).*(nxtarg.*nxtarg) + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*(nytarg.*nytarg))+...
               nu/(2*zk^2).*(hess(:, :, 1).*(tauxtarg.*tauxtarg) + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*(tauytarg.*tauytarg));
    
    submat = zeros(nkappa,2,nt,ns);
    submat(:,1,:,:) = reshape(firstbc,nkappa,1,nt,[]);
    submat(:,2,:,:) = reshape(secondbc,nkappa,1,nt,[]);
    submat = reshape(submat, [],ns);

case {'supported_plate_bcs_trx'}
    submat = zeros(nkappa,2,nt,4*ns);
   targnorm = targinfo.n;
   targtang = targinfo.d;

    nxtarg = repmat(reshape(targnorm(1,:),1,nt,1),nkappa,1,ns);
    nxtarg = reshape(nxtarg,[],ns);
    nytarg = repmat(reshape(targnorm(2,:),1,nt,1),nkappa,1,ns);
    nytarg = reshape(nytarg,[],ns);
    dx = repmat(reshape(targtang(1,:),1,nt,1),nkappa,1,ns);
    dx = reshape(dx,[],ns);
    dy = repmat(reshape(targtang(2,:),1,nt,1),nkappa,1,ns);
    dy = reshape(dy,[],ns);
    % dx = targinfo.d(1,:).';
    % dy = targinfo.d(2,:).';
    ds = sqrt(dx.*dx+dy.*dy);
    tauxtarg = (dx./ds);                                                                       % normalization
    tauytarg = (dy./ds);

    nu = varargin{1};
    
    [val, grad, hess,third, fourth, fifth] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,1);  
    
    firstbc = 1/(2*zk^2).*val ;
    
    secondbc = 1/(2*zk^2).*(hess(:, :, 1).*(nxtarg.*nxtarg) + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*(nytarg.*nytarg))+...
               nu/(2*zk^2).*(hess(:, :, 1).*(tauxtarg.*tauxtarg) + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*(tauytarg.*tauytarg));
    
    submat(:,1,:,1:4:end) = reshape(firstbc,nkappa,1,nt,[]);
    submat(:,2,:,1:4:end) = reshape(secondbc,nkappa,1,nt,[]);


    firstbc = 1/(2*zk^2).*grad(:, :, 1) ;
    
    secondbc = 1/(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg) + third(:, :, 2).*(2*nxtarg.*nytarg) + third(:, :, 3).*(nytarg.*nytarg))+...
               nu/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg) + third(:, :, 2).*(2*tauxtarg.*tauytarg) + third(:, :, 3).*(tauytarg.*tauytarg));
    
    submat(:,1,:,2:4:end) = -reshape(firstbc,nkappa,1,nt,[]);
    submat(:,2,:,2:4:end) = -reshape(secondbc,nkappa,1,nt,[]);

    firstbc = 1/(2*zk^2).*third(:, :, 1) ;
    
    secondbc = 1/(2*zk^2).*(fourth(:, :, 1).*(nxtarg.*nxtarg) + fourth(:, :, 2).*(2*nxtarg.*nytarg) + fourth(:, :, 3).*(nytarg.*nytarg))+...
               nu/(2*zk^2).*(fourth(:, :, 1).*(tauxtarg.*tauxtarg) + fourth(:, :, 2).*(2*tauxtarg.*tauytarg) + fourth(:, :, 3).*(tauytarg.*tauytarg));
    
    submat(:,1,:,3:4:end) = reshape(firstbc,nkappa,1,nt,[]);
    submat(:,2,:,3:4:end) = reshape(secondbc,nkappa,1,nt,[]);

        firstbc = 1/(2*zk^2).*fourth(:, :, 1) ;
    
    secondbc = 1/(2*zk^2).*(fifth(:, :, 1).*(nxtarg.*nxtarg) + fifth(:, :, 2).*(2*nxtarg.*nytarg) + fifth(:, :, 3).*(nytarg.*nytarg))+...
               nu/(2*zk^2).*(fifth(:, :, 1).*(tauxtarg.*tauxtarg) + fifth(:, :, 2).*(2*tauxtarg.*tauytarg) + fifth(:, :, 3).*(tauytarg.*tauytarg));
    
    submat(:,1,:,4:4:end) = -reshape(firstbc,nkappa,1,nt,[]);
    submat(:,2,:,4:4:end) = -reshape(secondbc,nkappa,1,nt,[]);

    submat = reshape(submat, [],4*ns);

% kernels for the supported plate integral equation, note there is no
% direct free-space analoge due to the kernel splitting
case {'supported_plate'}
    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    srcd2 = srcinfo.d2;
    targnorm = targinfo.n;
    targtang = targinfo.d;
    coefs = varargin{1};
    nu = coefs(1);

    if length(varargin) > 1
        nsub = varargin{2};
    else 
        nsub = 0;
    end

    nx = repmat(srcnorm(1,:),nkappa*nt,1);
    ny = repmat(srcnorm(2,:),nkappa*nt,1);
    
    nxtarg = repmat(reshape(targnorm(1,:),1,nt,1),nkappa,1,ns);
    nxtarg = reshape(nxtarg,[],ns);
    nytarg = repmat(reshape(targnorm(2,:),1,nt,1),nkappa,1,ns);
    nytarg = reshape(nytarg,[],ns);
        
    dx = repmat(srctang(1,:),nkappa*nt,1);
    dy = repmat(srctang(2,:),nkappa*nt,1);
    
    ds = sqrt(dx.*dx+dy.*dy);
    
    taux = dx./ds;
    tauy = dy./ds;

    dx1 = repmat(reshape(targtang(1,:),1,nt,1),nkappa,1,ns);
    dx1 = reshape(dx1,[],ns);
    dy1 = repmat(reshape(targtang(2,:),1,nt,1),nkappa,1,ns);
    dy1 = reshape(dy1,[],ns);

    
    ds1 = sqrt(dx1.*dx1+dy1.*dy1); 
    
    tauxtarg = dx1./ds1;
    tauytarg = dy1./ds1;
    
    dx1 = repmat(srctang(1,:),nkappa*nt,1);
    dy1 = repmat(srctang(2,:),nkappa*nt,1);
    
    d2x1 = repmat(srcd2(1,:),nkappa*nt,1);
    d2y1 = repmat(srcd2(2,:),nkappa*nt,1);
    
    denom = sqrt(dx1.^2+dy1.^2).^3;
    numer = dx1.*d2y1-d2x1.*dy1;
    
    curv = numer./denom;   
    
    kp = repmat(srcinfo.data(1,:),nkappa*nt,1);
    
    a1 = 2-nu;
    a2 = (-1+nu)*(7+nu)/(3 - nu);
    a3 = (1-nu)*(3+nu)/(1+nu);
    
    [~, grad, hess, third, fourth,fifth] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,0,nsub);  
    
    K11 = -1/(2*zk^2)*(third(:,:,1).*nx.^3 + 3*third(:,:,2).*nx.^2.*ny + 3*third(:,:,3).*nx.*ny.^2 + third(:,:,4).*ny.^3) + ...
         -a1/(2*zk^2)*(third(:,:,1).*nx.*taux.^2 + third(:,:,2).*(ny.*taux.^2 + 2*nx.*taux.*tauy) + third(:,:,3).*(nx.*tauy.^2 + 2*ny.*taux.*tauy) + third(:,:,4).*ny.*tauy.^2) + ...
         a2*curv./(2*zk^2).*(hess(:,:,1).*nx.^2 + 2*hess(:,:,2).*nx.*ny + hess(:,:,3).*ny.^2) + ...
         -a3*kp./(2*zk^2).*(grad(:,:,1).*taux + grad(:,:,2).*tauy);
    
    K12 = -1/(2*zk^2).*(grad(:, :, 1).*nx + grad(:, :, 2).*ny);
           
    K22 = -1/(2*zk^2).*(third(:,:,1).*(nx.*nxtarg.^2) + third(:,:,2).*(2*nx.*nxtarg.*nytarg + ny.*nxtarg.^2) + ...
         + third(:,:,3).*(nx.*nytarg.^2 + 2*ny.*nxtarg.*nytarg) + third(:,:,4).*(ny.*nytarg.^2))+...
         -nu./(2*zk^2).*(third(:,:,1).*(nx.*tauxtarg.^2) + third(:,:,2).*(2*nx.*tauxtarg.*tauytarg + ny.*tauxtarg.^2) + ...
         + third(:,:,3).*(nx.*tauytarg.^2 + 2*ny.*tauxtarg.*tauytarg) + third(:,:,4).*(ny.*tauytarg.^2)); 
    
    
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
         a2.*curv./(2*zk^2).*(fourth(:, :, 1).*(nxtarg.*nxtarg.*nx.*nx) + ...
          fourth(:, :, 2).*(2*nx.^2.*nxtarg.*nytarg + 2*nx.*nxtarg.^2.*ny) + ...
          fourth(:, :, 3).*(nxtarg.^2.*ny.^2 + 4*nx.*nxtarg.*ny.*nytarg + nx.^2.*nytarg.^2) +...
          fourth(:, :, 4).*(2*ny.*nytarg.*(nxtarg.*ny + nx.*nytarg))+...
          fourth(:, :, 5).*(nytarg.*nytarg.*ny.*ny)) + ... % kappa*G_{nx nx ny ny}
         nu*a2.*curv./(2*zk^2).*(fourth(:, :, 1).*(tauxtarg.*tauxtarg.*nx.*nx) + ...
          fourth(:, :, 2).*(2*nx.^2.*tauxtarg.*tauytarg + 2*nx.*tauxtarg.^2.*ny) + ...
          fourth(:, :, 3).*(tauxtarg.^2.*ny.^2 + 4*nx.*tauxtarg.*ny.*tauytarg + nx.^2.*tauytarg.^2) +...
          fourth(:, :, 4).*(2*ny.*tauytarg.*(tauxtarg.*ny + nx.*tauytarg))+...
          fourth(:, :, 5).*(tauytarg.*tauytarg.*ny.*ny)) + ... % kappa*nu*G_{nx nx ny ny}
        -a3.*kp./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
            third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
            third(:, :, 4).*(nytarg.*nytarg.*tauy)) + ... % kp*G_{nx nx tauy}
         -nu*a3.*kp./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
            third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
            third(:, :, 4).*(tauytarg.*tauytarg.*tauy)) ;  % kp*nu*G_{taux taux tauy}
    
    submat = zeros(nkappa,2,nt,2*ns);
    submat(:,1,:,1:2:2*ns) = reshape(K11,nkappa,1,nt,[]);
    submat(:,1,:,2:2:2*ns) = reshape(K12,nkappa,1,nt,[]);
    submat(:,2,:,1:2:2*ns) = reshape(K21,nkappa,1,nt,[]);
    submat(:,2,:,2:2:2*ns) = reshape(K22,nkappa,1,nt,[]);
    submat = reshape(submat, [],2*ns);


    
% supported plate kernels for plotting
case {'supported_plate_eval'}
    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    srcd2 = srcinfo.d2;
    coefs = varargin{1};
    nu = coefs(1);
    % nsub = varargin{2};
    
    nx = repmat(srcnorm(1,:),nkappa*nt,1);
    ny = repmat(srcnorm(2,:),nkappa*nt,1);

   dx = repmat(srctang(1,:),nkappa*nt,1);
   dy = repmat(srctang(2,:),nkappa*nt,1);
    
    ds = sqrt(dx.*dx+dy.*dy);
    
    taux = dx./ds;
    tauy = dy./ds;
    
    dx1 = repmat(srctang(1,:),nkappa*nt,1);
    dy1 = repmat(srctang(2,:),nkappa*nt,1);
    
    d2x1 = repmat(srcd2(1,:),nkappa*nt,1);
    d2y1 = repmat(srcd2(2,:),nkappa*nt,1);
    
    denom = sqrt(dx1.^2+dy1.^2).^3;
    numer = dx1.*d2y1-d2x1.*dy1;
    
    curv = numer./denom; 
    
    kp = repmat(srcinfo.data(1,:),nkappa*nt,1);
    
    a1 = 2-nu;
    a2 = (-1+nu)*(7+nu)/(3 - nu);
    a3 = (1-nu)*(3+nu)/(1+nu);
    
    [~, grad, hess, third] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,0);  
    % [val,grad] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,0);  
    
    K1 = -1/(2*zk^2)*(third(:,:,1).*nx.^3 + 3*third(:,:,2).*nx.^2.*ny + 3*third(:,:,3).*nx.*ny.^2 + third(:,:,4).*ny.^3) + ...
         -a1/(2*zk^2)*(third(:,:,1).*nx.*taux.^2 + third(:,:,2).*(ny.*taux.^2 + 2*nx.*taux.*tauy) + third(:,:,3).*(nx.*tauy.^2 + 2*ny.*taux.*tauy) + third(:,:,4).*ny.*tauy.^2) + ...
         a2*curv./(2*zk^2).*(hess(:,:,1).*nx.^2 + 2*hess(:,:,2).*nx.*ny + hess(:,:,3).*ny.^2) + ...
         -a3*kp./(2*zk^2).*(grad(:,:,1).*taux + grad(:,:,2).*tauy);

    K2 = -1/(2*zk^2).*(grad(:,:,1).*nx + grad(:,:,2).*ny);

    submat = zeros(nkappa*nt,2*ns);

    submat(:,1:2:end) = K1;
    submat(:,2:2:end) = K2;


case {'supported_plate_eval_trx'}

    submat = zeros(nkappa,4,nt,2*ns);
    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    srcd2 = srcinfo.d2;
    coefs = varargin{1};
    nu = coefs(1);
    
    nx = repmat(srcnorm(1,:),nkappa*nt,1);
    ny = repmat(srcnorm(2,:),nkappa*nt,1);

   dx = repmat(srctang(1,:),nkappa*nt,1);
   dy = repmat(srctang(2,:),nkappa*nt,1);
    
    ds = sqrt(dx.*dx+dy.*dy);
    
    taux = dx./ds;
    tauy = dy./ds;
    
    dx1 = repmat(srctang(1,:),nkappa*nt,1);
    dy1 = repmat(srctang(2,:),nkappa*nt,1);
    
    d2x1 = repmat(srcd2(1,:),nkappa*nt,1);
    d2y1 = repmat(srcd2(2,:),nkappa*nt,1);
    
    denom = sqrt(dx1.^2+dy1.^2).^3;
    numer = dx1.*d2y1-d2x1.*dy1;
    
    curv = numer./denom; 
    
    kp = repmat(srcinfo.data(1,:),nkappa*nt,1);
    
    a1 = 2-nu;
    a2 = (-1+nu)*(7+nu)/(3 - nu);
    a3 = (1-nu)*(3+nu)/(1+nu);
    
    [~, grad, hess, third,fourth,fifth,sixth] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,1);  
    
    K1 = -1/(2*zk^2)*(third(:,:,1).*nx.^3 + 3*third(:,:,2).*nx.^2.*ny + 3*third(:,:,3).*nx.*ny.^2 + third(:,:,4).*ny.^3) + ...
         -a1/(2*zk^2)*(third(:,:,1).*nx.*taux.^2 + third(:,:,2).*(ny.*taux.^2 + 2*nx.*taux.*tauy) + third(:,:,3).*(nx.*tauy.^2 + 2*ny.*taux.*tauy) + third(:,:,4).*ny.*tauy.^2) + ...
         a2*curv./(2*zk^2).*(hess(:,:,1).*nx.^2 + 2*hess(:,:,2).*nx.*ny + hess(:,:,3).*ny.^2) + ...
         -a3*kp./(2*zk^2).*(grad(:,:,1).*taux + grad(:,:,2).*tauy);

    K2 = -1/(2*zk^2).*(grad(:,:,1).*nx + grad(:,:,2).*ny);

    submat(:,1,:,1:2:end) = reshape(K1,nkappa,1,nt,[]);
    submat(:,1,:,2:2:end) = reshape(K2,nkappa,1,nt,[]);

    K1 = -1/(2*zk^2)*(fourth(:,:,1).*nx.^3 + 3*fourth(:,:,2).*nx.^2.*ny + 3*fourth(:,:,3).*nx.*ny.^2 + fourth(:,:,4).*ny.^3) + ...
         -a1/(2*zk^2)*(fourth(:,:,1).*nx.*taux.^2 + fourth(:,:,2).*(ny.*taux.^2 + 2*nx.*taux.*tauy) + fourth(:,:,3).*(nx.*tauy.^2 + 2*ny.*taux.*tauy) + fourth(:,:,4).*ny.*tauy.^2) + ...
         a2*curv./(2*zk^2).*(third(:,:,1).*nx.^2 + 2*third(:,:,2).*nx.*ny + third(:,:,3).*ny.^2) + ...
         -a3*kp./(2*zk^2).*(hess(:,:,1).*taux + hess(:,:,2).*tauy);

    K2 = -1/(2*zk^2).*(hess(:,:,1).*nx + hess(:,:,2).*ny);

    submat(:,2,:,1:2:end) = reshape(K1,nkappa,1,nt,[]);
    submat(:,2,:,2:2:end) = reshape(K2,nkappa,1,nt,[]);

    K1 = -1/(2*zk^2)*(fifth(:,:,1).*nx.^3 + 3*fifth(:,:,2).*nx.^2.*ny + 3*fifth(:,:,3).*nx.*ny.^2 + fifth(:,:,4).*ny.^3) + ...
         -a1/(2*zk^2)*(fifth(:,:,1).*nx.*taux.^2 + fifth(:,:,2).*(ny.*taux.^2 + 2*nx.*taux.*tauy) + fifth(:,:,3).*(nx.*tauy.^2 + 2*ny.*taux.*tauy) + fifth(:,:,4).*ny.*tauy.^2) + ...
         a2*curv./(2*zk^2).*(fourth(:,:,1).*nx.^2 + 2*fourth(:,:,2).*nx.*ny + fourth(:,:,3).*ny.^2) + ...
         -a3*kp./(2*zk^2).*(third(:,:,1).*taux + third(:,:,2).*tauy);

    K2 = -1/(2*zk^2).*(third(:,:,1).*nx + third(:,:,2).*ny);

    submat(:,3,:,1:2:end) = reshape(K1,nkappa,1,nt,[]);
    submat(:,3,:,2:2:end) = reshape(K2,nkappa,1,nt,[]);

    K1 = -1/(2*zk^2)*(sixth(:,:,1).*nx.^3 + 3*sixth(:,:,2).*nx.^2.*ny + 3*sixth(:,:,3).*nx.*ny.^2 + sixth(:,:,4).*ny.^3) + ...
         -a1/(2*zk^2)*(sixth(:,:,1).*nx.*taux.^2 + sixth(:,:,2).*(ny.*taux.^2 + 2*nx.*taux.*tauy) + sixth(:,:,3).*(nx.*tauy.^2 + 2*ny.*taux.*tauy) + sixth(:,:,4).*ny.*tauy.^2) + ...
         a2*curv./(2*zk^2).*(fifth(:,:,1).*nx.^2 + 2*fifth(:,:,2).*nx.*ny + fifth(:,:,3).*ny.^2) + ...
         -a3*kp./(2*zk^2).*(third(:,:,1).*taux + third(:,:,2).*tauy);

    K2 = -1/(2*zk^2).*(third(:,:,1).*nx + third(:,:,2).*ny);

    submat(:,4,:,1:2:end) = reshape(K1,nkappa,1,nt,[]);
    submat(:,4,:,2:2:end) = reshape(K2,nkappa,1,nt,[]);

    submat = reshape(submat, [], 2*ns);

end

if (ising == 1) && ~contains(type,'trx')
    ishape = size(submat);
    submat = reshape(submat,nkappa,ishape(1)/nkappa,ishape(2));
    submat = submat + reshape(chnk.flex2d.kern(zk,srcinfo,targinfo,type,varargin{:}),1,ishape(1)/nkappa,ishape(2));
    submat = reshape(submat,ishape);
end
end


