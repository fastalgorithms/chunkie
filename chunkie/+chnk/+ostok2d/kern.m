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
    [~,grad] = chnk.lap2d.green(k,src,targ);
  for i = 1:2
      for j = 1:2
           for l = 1:2
              d_ij = 1-abs(i-j);   %%%%%%%% kronecker delta for i and j
              d_il = 1-abs(i-l);
              d_lj = 1-abs(l-j);
              T(i,j,l) = -d_il*grad(:,:,j) - d_ij*(der3(:,:,2+l-2) + der3(:,:,4+l-2))...
                  -d_lj*(der3(:,:,2+i-2) + der3(:,:,4+i-2)) + 2*der3(:,:,i+j+l-2);
          end
      end
  end
  submat = T;
end
