
function submat = glapkern(src,targ,srctau,targtau,type)

[~,ns] = size(src);
[~,nt] = size(targ);

if strcmpi(type,'d')
    [~,grad] = glapfun(src,targ);
    nx = repmat(srctau(2,:),nt,1);
    ny = repmat(-srctau(1,:),nt,1);
    submat = -(grad(:,:,1).*nx + grad(:,:,2).*ny);
end

if strcmpi(type,'sprime')
    [~,grad] = glapfun(src,targ);
    nx = repmat((targtau(2,:)).',1,ns);
    ny = repmat(-(targtau(1,:)).',1,ns);

    submat = (grad(:,:,1).*nx + grad(:,:,2).*ny);
end

if strcmpi(type,'s')
    submat = glapfun(src,targ);
end