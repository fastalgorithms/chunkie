
function submat = kern(zk,src,targ,srctau,targtau,type,varargin)
%CHNK.MFIE.KERN F.T. in z of MFIE kernel
% this function evaluates the matrix corresponding to the kernel
% for the given sources and targets
% 

[~,ns] = size(src);
[~,nt] = size(targ);

if strcmpi(type,'xx')
    xi = varargin{1};
    [~,grad] = chnk.helm2d.green(zk,src,targ);
%     nx = repmat(srctau(2,:),nt,1);
%     ny = repmat(-srctau(1,:),nt,1);
    submat = grad(:,:,1).*srctau(2,:) - grad(:,:,2).*srctau(1,:);
end

if strcmpi(type,'zz')
    [~,grad] = chnk.helm2d.green(zk,src,targ);
    nx = repmat((targtau(2,:)).',1,ns);
    ny = repmat(-(targtau(1,:)).',1,ns);

    submat = -(grad(:,:,1).*nx + grad(:,:,2).*ny);
end
