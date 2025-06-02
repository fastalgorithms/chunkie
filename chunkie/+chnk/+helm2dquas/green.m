function [val,grad,hess] = green(src,targ,zk,kappa,d,sn,l)
%CHNK.HELM2DQUAS.GREEN evaluate the quasiperiodic helmholtz green's function
% for the given sources and targets
%
% Input:
%   zk - wavenumber
%   kappa - quasiperiodic parameters
%   d - period
%   sn - precomputed lattice sum integrals 
%       (see chnk.helm2dquas.latticecoefs)
%
% see also CHNK.HELM2DQUAS.KERN
[~,nsrc] = size(src);
[~,ntarg] = size(targ);


xs = repmat(src(1,:),ntarg,1);
ys = repmat(src(2,:),ntarg,1);

xt = repmat(targ(1,:).',1,nsrc);
yt = repmat(targ(2,:).',1,nsrc);

rx = xt-xs;
ry = yt-ys;

nx = fix(rx/d);
rx = rx - nx*d;


rx2 = rx.*rx;
ry2 = ry.*ry;

r2 = rx2+ry2;

r = sqrt(r2);

rx = rx(:);
ry = ry(:);
r = r(:);

npt = size(r,1);

ythresh = d/2;
% ythresh = d/10;
iclose = abs(ry) < ythresh;
ifar = ~iclose;

rxfar = rx(ifar);
ryfar = ry(ifar);

rxclose = rx(iclose);
ryclose = ry(iclose);
rclose = r(iclose);

nkappa = length(kappa);

val = zeros(nkappa,npt,1);
if nargout > 1
grad = zeros(nkappa,npt,2);
end
if nargout > 2
hess = zeros(nkappa,npt,3);
end


tol = 1e-10;
Lbd = sqrt((log(tol))^2/real(ythresh)^2 + real(zk)^2);

rxfar = rxfar.';
ryfar = ryfar.';
if ~isempty(ryfar)
M = ceil(Lbd*d/(2*pi));
ms = reshape((-M:M),1,1,[]);
xi_m = kappa(:) + 2*pi/d*ms;

% beta = sqrt((xi_m.^2-zk^2));
beta = sqrt(1i*(xi_m-zk)).*sqrt(-1i*(xi_m+zk));

% fhat = exp(-beta.*sqrt(ryfar.^2))./beta.*exp(1i*xi_m.*rxfar)/2;
fhat = exp(-beta.*sqrt(ryfar.^2) + 1i*xi_m.*rxfar)./(2*beta);
val(:,ifar,:) = sum(fhat,3)/(d);
if nargout > 1
grad(:,ifar,1) = sum(1i*xi_m.*fhat,3)/d;
grad(:,ifar,2) = sum(-beta.*(sqrt(ryfar.^2)./ryfar).*fhat,3)/d;
% grad(:,ifar,:) =[gx,gy];
end

if nargout >2
hess(:,ifar,1) = sum(-xi_m.^2.*fhat,3)/d;
hess(:,ifar,2) = sum(-1i*xi_m.*beta.*(sqrt(ryfar.^2)./ryfar).*fhat,3)/d;
hess(:,ifar,3) = sum((beta.*(sqrt(ryfar.^2)./ryfar)).^2.*fhat,3)/d;
% hess(:,ifar,:) = [hxx,hxy,hyy];
end


end

alpha = (exp(1i*kappa(:)*d));

val_near=0;
grad_near = zeros(1,1,2);
hess_near = zeros(1,1,3);
if ~isempty(rxclose)
    for i = -l:l
        rxi = rxclose - i*d;
        if nargout>2
        [vali,gradi,hessi] = chnk.helm2d.green(zk,[0;0],[rxi.';ryclose.']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        hessi = reshape(hessi,1,[],3);
        val_near = val_near + vali.*alpha.^i;
        grad_near = grad_near + gradi.*alpha.^i;
        hess_near = hess_near + hessi.*alpha.^i;
        elseif nargout > 1
        [vali,gradi] = chnk.helm2d.green(zk,[0;0],[rxi.';ryclose.']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        val_near = val_near + vali.*alpha.^i;
        grad_near = grad_near + gradi.*alpha.^i;
        else
        vali = chnk.helm2d.green(zk,[0;0],[rxi.';ryclose.']);
        vali = reshape(vali,1,[],1);
        val_near = val_near + vali.*alpha.^i;
        end
    end
    
    % sn = sn.';
    N = size(sn,2)-1;
    ns = (0:N);
    ns_use = (0:N+2);
    Js = zeros(length(rclose),N+3);

    if length(rclose) < N+3
        for i = 1:length(rclose)
        Js(i,:) = besselj(ns_use,zk*rclose(i));
        end
    else
        for i = 1:length(ns_use)
        Js(:,i) = besselj(ns_use(i),zk*rclose);
        end
    end
    eip = (rxclose+1i*ryclose)./rclose;
    eipn = reshape(eip.^ns,1,[], N+1);
    cs = (eipn+1./eipn)/2;
    
    Js = reshape(Js,1,[],N+3);
    sn = reshape(sn,nkappa, 1, N+1);
    
    val_far = 0.25*1i*Js(:,:,1).*sn(:,:,1) + 0.5*1i*sum(sn(:,:,2:end).*Js(:,:,2:end-2).*cs(:,:,2:end),3);
    val(:,iclose) = val_near+val_far;
    
    
    
    
    if nargout >1
        DJs = cat(3,-Js(:,:,2),.5*(Js(:,:,1:end-3)-Js(:,:,3:end-1)))*zk;
        ss = (eipn-1./eipn)/2i;
            
        grad_far_p = 0.25*1i*DJs(:,:,1).*sn(:,:,1) + 0.5*1i*sum(sn(:,:,2:end).*DJs(:,:,2:end).*cs(:,:,2:end),3);
        grad_far_t = (0.5*1i*sum(-reshape((1:N),1,1,[]).*sn(:,:,2:end).*Js(:,:,2:end-2).*ss(:,:,2:end),3))./rclose.';
        
        grad_far = cat(3,cs(:,:,2).*grad_far_p - ss(:,:,2).*grad_far_t, ss(:,:,2).*grad_far_p + cs(:,:,2).*grad_far_t);
        grad(:,iclose,:) = grad_near + grad_far; 
    end
    if nargout > 2
        DDJs = cat(3,.5*(Js(:,:,3)-Js(:,:,1)),.25*(Js(:,:,4)-3*Js(:,:,2)),.25*(Js(:,:,1:end-4)-2*Js(:,:,3:end-2)+Js(:,:,5:end)))*zk^2;
        rclose = rclose.';
        rxclose = rxclose.';
        ryclose = ryclose.';
        ns = reshape(ns,1,1,[]);

        tmp_n = rclose.^(-4).*(-ns.*ryclose.*Js(:,:,1:end-2).*(ns.*ryclose.*cs+2*rxclose.*ss)+ ...
            rclose.*ryclose.*(ryclose.*cs + 2*ns.*rxclose.*ss).*DJs+ ...
            rclose.^2.*rxclose.^2.*cs.*DDJs);
        hess_far_xx = 0.25*1i*tmp_n(:,:,1).*sn(:,:,1)+.5*1i*sum(sn(:,:,2:end).*tmp_n(:,:,2:end),3);

        tmp_n = ns.*Js(:,:,1:end-2).*(ns.*rxclose.*ryclose.*cs+(rxclose.^2-ryclose.^2).*ss).*rclose.^(-4) ...
            +rclose.^(-3).*(-(rxclose.*ryclose.*cs+ns.*(rxclose.^2-ryclose.^2).*ss).*DJs+...
            rclose.*rxclose.*ryclose.*cs.*DDJs);


        hess_far_xy = 0.25*1i*tmp_n(:,:,1).*sn(:,:,1)+.5*1i*sum(sn(:,:,2:end).*tmp_n(:,:,2:end),3);

        tmp_n = rclose.^(-4).*(-ns.*rxclose.*Js(:,:,1:end-2).*(ns.*rxclose.*cs-2*ryclose.*ss)+ ...
            rclose.*rxclose.*(rxclose.*cs - 2*ns.*ryclose.*ss).*DJs+ ...
            rclose.^2.*ryclose.^2.*cs.*DDJs);
        hess_far_yy = 0.25*1i*tmp_n(:,:,1).*sn(:,:,1)+.5*1i*sum(sn(:,:,2:end).*tmp_n(:,:,2:end),3);

        hess_far = cat(3,hess_far_xx, hess_far_xy, hess_far_yy);


        hess(:,iclose,:) = hess_near + hess_far;
    end
end

quasi_phase = exp(1i*kappa(:)*nx(:).'*d);

val = reshape(quasi_phase.*val,nkappa*ntarg,nsrc);
if nargout>1
grad = reshape(quasi_phase.*grad,nkappa*ntarg,nsrc,2);
end
if nargout>2
hess = reshape(quasi_phase.*hess,nkappa*ntarg,nsrc,3);
end
% end

end

