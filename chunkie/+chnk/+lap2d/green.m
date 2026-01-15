%%%%%%   LAPLACE GREEN   %%%%%%%%
function [val,grad,hess] = green(src,targ,nolog)

if (nargin < 3)
    nolog = false;
end

[~,ns] = size(src);
[~,nt] = size(targ);

rx = bsxfun(@minus,targ(1,:).',src(1,:));
ry = bsxfun(@minus,targ(2,:).',src(2,:));
r2 = rx.^2+ry.^2;


if (nargout > 0) 
    
    if nolog
        val = [];
    else
        val = -log(r2)/(4.0*pi);
    end
    
end

if nargout > 1 

    grad = zeros(nt,ns,2);

    grad(:,:,1) = -rx./(2.0*pi*r2);
    grad(:,:,2) = -ry./(2.0*pi*r2);
    
end

if nargout > 2 
    rx2 = rx.*rx;
    ry2 = ry.*ry;


    r4 = r2.*r2;

    hess = zeros(nt,ns,3);

    hess(:,:,1) = rx2./(pi*r4)-1.0./(2.0*pi*r2);
    hess(:,:,2) = rx.*ry./(pi*r4);
    hess(:,:,3) = ry2./(pi*r4)-1.0./(2.0*pi*r2);
end
