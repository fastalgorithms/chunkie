%%%%%%%%%   KERNEL FOR LAYER POTENTIALS FOR OSCILLATORY STOKES %%%%%%
function varargout = kern(k,srcinfo,targinfo,type,varargin)
% for type 's' submat gives the usual GREEN's function for oscillatory
% STOKES operator 
% []
% for type 'd' submat gives the Double Layer for oscillatory
% STOKES operator 
% []
% for type 'c' submat gives the Combined field for oscillatory
% STOKES operator


src = srcinfo.r(:,:);
targ = targinfo.r(:,:);
[~, nt] = size(targ);
[~, ns] = size(src);

% Single Layer 

if strcmpi(type,'s')
  
  [~,~, hess] = chnk.flex2d.helmdiffgreen(k,src,targ);
  hess = hess/(k*k);
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
  
[~,~,~,der3] = chnk.flex2d.helmdiffgreen(k,src,targ);
der3 = der3/(k*k);
[~, nt] = size(targinfo.r(:,:));
[~, ns] = size(srcinfo.r(:,:));


nx = repmat(srcinfo.n(1,:),nt,1);
ny = repmat(srcinfo.n(2,:),nt,1);

[~, grad] = chnk.lap2d.green(src, targ);

Kxx = -((grad(:,:,1) - 2*der3(:,:,3)).*nx ...
         + (der3(:,:,2) - der3(:,:,4)).*ny);

Kxy = -((grad(:,:,2) + 2*der3(:,:,2)).*nx ...
          + (der3(:,:,3) - der3(:,:,1)).*ny);

Kyx = -((der3(:,:,2) - der3(:,:,4)).*nx ...
          + (grad(:,:,1) + 2*der3(:,:,3)).*ny);

Kyy = -((der3(:,:,3) - der3(:,:,1)).*nx ...
          + (grad(:,:,2) - 2*der3(:,:,2)).*ny);


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

if strcmpi(type,'c')
        coefs = varargin{1};
        [Sxx, Sxy, Syy] = chnk.ostok2d.kern(k, srcinfo, targinfo, 's');
        [Dxx, Dxy, Dyx, Dyy] = chnk.ostok2d.kern(k, srcinfo,targinfo, 'd');
        Kxx = coefs(1)*Dxx + coefs(2)*Sxx;
        Kyy = coefs(1)*Dyy + coefs(2)*Syy;
        Kxy = coefs(1)*Dxy + coefs(2)*Sxy;
        Kyx = coefs(1)*Dyx + coefs(2)*Sxy;

        if ( nargout == 1 )
            % Interleave
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

end