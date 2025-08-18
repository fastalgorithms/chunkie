%%%%%%%%%   KERNEL FOR LAYER POTENTIALS FOR OSCILLATORY STOKES %%%%%%
function varargout = kern(k,srcinfo,targinfo,type,varargin)
% for type 's' submat gives the usual GREEN's function for oscillatory
% STOKES operator
%
% for type 'd' submat gives the [T_ijl]

  
src = srcinfo.r(:,:);
targ = targinfo.r(:,:);
[~, nt] = size(targ);
[~, ns] = size(src);




% single layer 

if strcmpi(type,'s')
  [~, ~, hess] = chnk.obihar2d.green(k, src, targ);
  Kxx = -hess(:,:,3);
  Kxy = hess(:,:,2);
  Kyy = -hess(:,:,1);
  if nargout == 1
      K = zeros(2*nt, 2*ns);
      K(1:2:end, 1:2:end) = Kxx;
      K(1:2:end, 2:2:end) = Kxy;
      K(2:2:end, 1:2:end) = Kxy;
      K(2:2:end, 2:2:end) = Kyy;
      varargout = {K};
  else
      varargout = {Kxx, Kxy, Kyy};
  end
end

% double layer 

if strcmpi(type,'d')
  
[~,~,~,der3] = chnk.obihar2d.green(k,src,targ);
[~, nt] = size(targinfo.r(:,:));
[~, ns] = size(srcinfo.r(:,:));


nx = repmat(srcinfo.n(1,:),nt,1);
ny = repmat(srcinfo.n(2,:),nt,1);

[~, grad] = chnk.lap2d.green(src, targ);

Kxx = -((grad(:,:,1) - 2*der3(:,:,3)).*nx + (der3(:,:,2) - der3(:,:,4)).*ny);

Kxy = -((grad(:,:,2) + 2*der3(:,:,2)).*nx + (der3(:,:,3) - der3(:,:,1)).*ny);

Kyx = -((der3(:,:,2) - der3(:,:,4)).*nx + (grad(:,:,1) + 2*der3(:,:,3)).*ny);

Kyy = -((der3(:,:,3) - der3(:,:,1)).*nx + (grad(:,:,2) - 2*der3(:,:,2)).*ny);



 if nargout == 1
      K = zeros(2*nt, 2*ns);     
      K(1:2:end, 1:2:end) = Kxx;
      K(1:2:end, 2:2:end) = Kyx;
      K(2:2:end, 1:2:end) = Kxy;
      K(2:2:end, 2:2:end) = Kyy;
      varargout = {K};
  else
      varargout = {Kxx, Kxy, Kyx, Kyy};
  end
 
end


% sprime layer 

if strcmpi(type,'sp')
  
[~,~,~,der3] = chnk.obihar2d.green(k,src,targ);
[~, nt] = size(targinfo.r(:,:));
[~, ns] = size(srcinfo.r(:,:));


nx = repmat(targinfo.n(1,:).',1,ns);
ny = repmat(targinfo.n(2,:).',1,ns);

[~, grad] = chnk.lap2d.green(src, targ);

Kxx = ((grad(:,:,1) - 2*der3(:,:,3)).*nx + (der3(:,:,2) - der3(:,:,4)).*ny);

Kxy = ((grad(:,:,2) + 2*der3(:,:,2)).*nx + (der3(:,:,3) - der3(:,:,1)).*ny);

Kyx = ((der3(:,:,2) - der3(:,:,4)).*nx + (grad(:,:,1) + 2*der3(:,:,3)).*ny);

Kyy = ((der3(:,:,3) - der3(:,:,1)).*nx + (grad(:,:,2) - 2*der3(:,:,2)).*ny);


 if nargout == 1
      K = zeros(2*nt, 2*ns);     
      K(1:2:end, 1:2:end) = Kxx;
      K(1:2:end, 2:2:end) = Kxy;
      K(2:2:end, 1:2:end) = Kyx;
      K(2:2:end, 2:2:end) = Kyy;
      varargout = {K};
  else
      varargout = {Kxx, Kxy, Kyx, Kyy};
  end
 
end

% sink velocity
% for sv, we have u = grad(phi) where phi = -1/(2*pi)*logr the green's function for laplacian
% and for the pressure we have p = delta - k^2/(2*pi)*logr   

if strcmpi(type, 'sinkv')
    l = k^2;
    [val, grad, hess] = chnk.lap2d.green(src, targ);

    %%%%%%%%% green's function = grad(Green_laplacian) %%%%%%

    % u1 = grad(:,:,1);
    % u2 = grad(:,:,2);
    % 
    % 
    % if nargout == 1
    %   K = zeros(2*nt, 2*ns);     
    %   K(1:2:end, 1:2:end) = u1;
    %   K(1:2:end, 2:2:end) = u2;
    %   K(2:2:end, 1:2:end) = u1;
    %   K(2:2:end, 2:2:end) = u2;
    %   varargout = {K};
    % else
    %   varargout = {u1, u2, u1, u2};
    % end

    %%%%%%%%%%%%%  sigma = -p.Id + (grad(u) + grad(u).')

    s11 = -l*val + hess(:,:,1);
    s12 = hess(:,:,2);
    s22 = -l*val + hess(:,:,3);
    % 
    % if nargout == 1
    %   K = zeros(2*nt, 2*ns);     
    %   K(1:2:end, 1:2:end) = s11;
    %   K(1:2:end, 2:2:end) = s12;
    %   K(2:2:end, 1:2:end) = s12;
    %   K(2:2:end, 2:2:end) = s22;
    %   varargout = {K};
    % else
    %   varargout = {s11, s12, s12, s22};
    % end


    %%%%%%%%%%%%%%%  traction = sigma.nu
    nx = repmat(srcinfo.n(1,:),nt,1);
    ny = repmat(srcinfo.n(2,:),nt,1);

    u1 = s11.*nx + s12.*ny;
    u2 = s22.*nx + s12.*ny;


    if nargout == 1
      K = zeros(2*nt, 2*ns);     
      K(1:2:end, 1:2:end) = u1;
      K(1:2:end, 2:2:end) = u2;
      K(2:2:end, 1:2:end) = u1;
      K(2:2:end, 2:2:end) = u2;
      varargout = {K};
    else
      varargout = {u1, u2, u1, u2};
    end

end