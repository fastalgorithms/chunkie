function [h] = newt_step(h,umesh,dmesh,qmesh,sig0,lam,opts)

level = 0.5;
step_fact = 1;
ww = qmesh.wts(:).';

if (nargin > 6)
if (isfield(opts,"levels"))
    level = opts.levels;
end

if (isfield(opts,"scales"))
    ww = ww.*opts.scales.';
end

if (isfield(opts,"step_fact"))
    step_fact = opts.step_fact;
end
end

dpx = dmesh.pseudo_normals(1,:).';
dpy = dmesh.pseudo_normals(2,:).';

rnx = qmesh.n(1,:).';
rny = qmesh.n(2,:).';

    rt = dmesh.r;
    rt(1,:) = rt(1,:) + (h.*dpx).';
    rt(2,:) = rt(2,:) + (h.*dpy).';

    [sig, sig_grad] = chnk.smoother.get_sigs(umesh, rt, sig0, lam);
    [val, grad, hess, hess_sig] = chnk.smoother.green(qmesh.r, rt, sig);
    gx = grad(:,:,1); 
    gy = grad(:,:,2);
    
    phi = -(gx*(rnx.*ww.') + gy*(rny.*ww.')) - level;
    
    
    h11 = hess(:,:,1,1).*ww;
    h12 = hess(:,:,1,2).*ww;
    h21 = hess(:,:,2,1).*ww;
    h22 = hess(:,:,2,2).*ww;
    
    h1sig = hess_sig(:,:,1).*ww;
    h2sig = hess_sig(:,:,2).*ww;

    dx1 = h11*rnx + h12*rny;
    dy1 = h21*rnx + h22*rny;
    
    dx2 = h1sig*rnx + h2sig*rny;
    dy2 = dx2;
    
    dsigx = sig_grad(:,1);
    dsigy = sig_grad(:,2);

    dx2 = dx2.*dsigx;
    dy2 = dy2.*dsigy;
    
    
    dphidh = (dx1 + dx2).*dpx + (dy1 + dy2).*dpy;
    h = h - step_fact*phi./dphidh;

end