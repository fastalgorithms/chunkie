function [val,grad,hess,hess_sig] = green(rs,rt,dsigt)

[~,ns] = size(rs);
[~,nt] = size(rt);

xs = repmat(rs(1,:),nt,1);
ys = repmat(rs(2,:),nt,1);

xt = repmat(rt(1,:).',1,ns);
yt = repmat(rt(2,:).',1,ns);
dst = repmat(dsigt.',1,ns);

dx = xt-xs;
dy = yt-ys;
dr = sqrt(dx.^2+dy.^2);

sz = size(dx);
dx = dx(:);
dy = dy(:);
dr = dr(:);
dst = dst(:);

da = sqrt(2.0d0)*dsigt;
z = dr./da;

iloc = find(z < 0.1);
ifar = find(z >= 0.1);

vs = zeros(ns*nt,1);
gs = zeros(ns*nt,2);
hs = zeros(ns*nt,2,2);
hsigs = zeros(ns*nt,2);

if (numel(iloc) ~=0)
    [vl,gl,hl,hsl] = chnk.smoother.gpsi_loc(dx,dy,dr,dst);
    vs(iloc) = vl;
    gs(iloc,:) = gl;
    hs(iloc,:,:) = hl;
    hsigs(iloc,:) = hsl;
end
if (numel(ifar) ~=0)
    [vf,gf,hf,hsf] = chnk.smoother.gpsi_all(dx,dy,dr,dsigt);
    vs(ifar) = vf;
    gs(ifar,:) = gf;
    hs(ifar,:,:) = hf;
    hsigs(ifar,:) = hsf;
end

val = reshape(vs,sz);
if (nargout > 1)
    grad = squeeze(reshape(gs,[sz,2]));
end
if (nargout > 2)
    hess = squeeze(reshape(hs,[sz,2,2]));
end
if (nargout > 3)
    hess_sig = squeeze(reshape(hsigs,[sz,2]));

end

end
