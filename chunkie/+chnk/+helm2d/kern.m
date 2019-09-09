
function submat = kern(zk,src,targ,srctau,targtau,type,varargin)
%CHNK.HELM2D.KERN standard Helmholtz kernels in 2D
% this function evaluates the matrix corresponding to the kernel
% for the given sources and targets
% 

[~,ns] = size(src);
[~,nt] = size(targ);

if strcmpi(type,'d')
    [~,grad] = chnk.helm2d.green(zk,src,targ);
    nx = repmat(srctau(2,:),nt,1);
    ny = repmat(-srctau(1,:),nt,1);
    submat = -(grad(:,:,1).*nx + grad(:,:,2).*ny);
end

if strcmpi(type,'sprime')
    [~,grad] = chnk.helm2d.green(zk,src,targ);
    nx = repmat((targtau(2,:)).',1,ns);
    ny = repmat(-(targtau(1,:)).',1,ns);

    submat = (grad(:,:,1).*nx + grad(:,:,2).*ny);
end

if strcmpi(type,'s')
    submat = chnk.helm2d.green(zk,src,targ);
end

if strcmpi(type,'c')
    eta = varargin{1};
    [s,grad] = chnk.helm2d.green(zk,src,targ);
    nx = repmat(srctau(2,:),nt,1);
    ny = repmat(-srctau(1,:),nt,1);
    submat = -(grad(:,:,1).*nx + grad(:,:,2).*ny) + 1i*eta*s;
end
