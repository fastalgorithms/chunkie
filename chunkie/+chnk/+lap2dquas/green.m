function [val,grad,hess] = green(src,targ,kappa,d,s0,sn,l,ising,nsub)
%CHNK.LAP2DQUAS.GREEN evaluate the quasiperiodic Helmholtz Green's function
% for the given sources and targets
%
% Input:
%   kappa - quasiperiodic parameters
%   d - period
%   s0 - precomputed lattice sum for n = 0 
%   sn - precomputed lattice sum for n >1 
%   l - number of periodic copies computed explicitly
%   ising - if set to 0, only include the periodic copies. If set to 1,
%       include the free-space part
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


rx = rx(:);
ry = ry(:);


nx = fix(rx/d);
rx = rx - nx*d;


rx2 = rx.*rx;
ry2 = ry.*ry;

r2 = rx2+ry2;

r = sqrt(r2);

npt = size(r,1);

ythresh = 1.5*d/2;
iclose = abs(ry) < ythresh;
ifar = ~iclose;

rxfar = rx(ifar);
ryfar = ry(ifar);

rxclose = rx(iclose);
ryclose = ry(iclose);
rclose = r(iclose);

nxclose = nx(iclose);
nptclose = size(rxclose, 1);

nkappa = length(kappa);

val = zeros(nkappa,npt,1);
if nargout > 1
grad = zeros(nkappa,npt,2);
end
if nargout > 2
hess = zeros(nkappa,npt,3);
end

if nargin < 9
    nsub = 0;
elseif nsub > l
    error('trying to subtract off too many copies')
end

zk = 0;
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

fhat = exp(-beta.*sqrt(ryfar.^2) + 1i*xi_m.*rxfar)./(2*beta);
val(:,ifar,:) = sum(fhat,3)/(d);
if nargout > 1
grad(:,ifar,1) = sum(1i*xi_m.*fhat,3)/d;
grad(:,ifar,2) = sum(-beta.*(sqrt(ryfar.^2)./ryfar).*fhat,3)/d;
end

if nargout >2
hess(:,ifar,1) = sum(-xi_m.^2.*fhat,3)/d;
hess(:,ifar,2) = sum(-1i*xi_m.*beta.*(sqrt(ryfar.^2)./ryfar).*fhat,3)/d;
hess(:,ifar,3) = sum((beta.*(sqrt(ryfar.^2)./ryfar)).^2.*fhat,3)/d;
end

end

alpha = (exp(1i*kappa(:)*d));

val_near= zeros(nkappa,nptclose);
grad_near = zeros(nkappa,nptclose,2);
hess_near = zeros(nkappa,nptclose,3);
ls = -l:l;
if ~isempty(rxclose)
    for i = ls
        if ising == 1
            iuse = true(nptclose,1);
        else
            iuse = ~ismember(nxclose, -i-nsub:-i+nsub); 
        end

        rxi = rxclose - i*d;
        if nargout>2
        [vali,gradi,hessi] = chnk.lap2d.green([0;0],[rxi.';ryclose.']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        hessi = reshape(hessi,1,[],3);
        val_near(:,iuse) = val_near(:,iuse) + vali(:,iuse).*alpha.^i;
        grad_near(:,iuse,:) = grad_near(:,iuse,:) + gradi(:,iuse,:).*alpha.^i;
        hess_near(:,iuse,:) = hess_near(:,iuse,:) + hessi(:,iuse,:).*alpha.^i;
        elseif nargout > 1
        [vali,gradi] = chnk.lap2d.green([0;0],[rxi.';ryclose.']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        val_near(:,iuse) = val_near(:,iuse) + vali(:,iuse).*alpha.^i;
        grad_near(:,iuse,:) = grad_near(:,iuse,:) + gradi(:,iuse,:).*alpha.^i;
        else
        vali = chnk.lap2d.green([0;0],[rxi.';ryclose.']);
        vali = reshape(vali,1,[],1);
        val_near(:,iuse) = val_near(:,iuse) + vali(:,iuse).*alpha.^i;
        end
    end

    N = size(sn,2);
    ns = 1:N;
    
    Rns = rclose.^ns;
    
    eip = (rxclose+1i*ryclose)./rclose;
    eipn = reshape(eip.^ns,1,[], N);
    cs = (eipn+1./eipn)/2;
    
    Rns = reshape(Rns,1,[],N);
    sn = reshape(sn,nkappa, 1, N);
    
    tmp = reshape(Rns.*cs,[],N);
    val_far = sn(:,:)*tmp.' + s0;
    
    ndiag = sum(rclose < 1e-14);
    if ndiag > 0
    val_far(:,rclose < 1e-14) = repmat(s0,1,ndiag); % diagonal replacement
    end

    val(:,iclose) = val_near+val_far;
    
    if nargout >1
        DRns = ns.*rclose.^(ns-1);
        DRns = reshape(DRns,1,[],N);
        ss = (eipn-1./eipn)/2i;
            
        tmp = reshape(DRns.*cs,[],N);
        grad_far_p = sn(:,:)*tmp.';
        tmp = reshape(Rns.*ss,[],N)./rclose;
        grad_far_t = ((-reshape((1:N),1,[]).*sn(:,:))*tmp.');
        
        grad_far = cat(3,cs(:,:,1).*grad_far_p - ss(:,:,1).*grad_far_t, ss(:,:,1).*grad_far_p + cs(:,:,1).*grad_far_t);
        
        if ndiag > 0
            grad_far(:,rclose < 1e-14,1) = repmat(sn(:,1),1,ndiag);
            grad_far(:,rclose < 1e-14,2) = 0;            
        end
        
        grad(:,iclose,:) = grad_near + grad_far; 
    end
    if nargout > 2
        DDRns = ns.*(ns-1).*rclose.^(ns-2);
        DDRns = reshape(DDRns,1,[],N);
        rclose = rclose.';
        rxclose = rxclose.';
        ryclose = ryclose.';
        ns = reshape(ns,1,1,[]);

        tmp_n = rclose.^(-4).*(-ns.*ryclose.*Rns(:,:,:).*(ns.*ryclose.*cs+2*rxclose.*ss)+ ...
            rclose.*ryclose.*(ryclose.*cs + 2*ns.*rxclose.*ss).*DRns+ ...
            rclose.^2.*rxclose.^2.*cs.*DDRns);
        tmp_n = reshape(tmp_n,[],N);
        hess_far_xx = sn(:,:)*tmp_n(:,:).';

        tmp_n = ns.*Rns(:,:,:).*(ns.*rxclose.*ryclose.*cs+(rxclose.^2-ryclose.^2).*ss).*rclose.^(-4) ...
            +rclose.^(-3).*(-(rxclose.*ryclose.*cs+ns.*(rxclose.^2-ryclose.^2).*ss).*DRns+...
            rclose.*rxclose.*ryclose.*cs.*DDRns);

        tmp_n = reshape(tmp_n,[],N);
        hess_far_xy = sn(:,:)*tmp_n(:,:).';

        tmp_n = rclose.^(-4).*(-ns.*rxclose.*Rns(:,:,:).*(ns.*rxclose.*cs-2*ryclose.*ss)+ ...
            rclose.*rxclose.*(rxclose.*cs - 2*ns.*ryclose.*ss).*DRns+ ...
            rclose.^2.*ryclose.^2.*cs.*DDRns);
        tmp_n = reshape(tmp_n,[],N);
        hess_far_yy = sn(:,:)*tmp_n(:,:).';

        hess_far = cat(3,hess_far_xx, hess_far_xy, hess_far_yy);

        if ndiag > 0
            hess_far(:,rclose < 1e-14,1) = repmat(2*sn(:,2),1,ndiag);
            hess_far(:,rclose < 1e-14,2) = 0; 
            hess_far(:,rclose < 1e-14,3) = repmat(-2*sn(:,2),1,ndiag);            
        end

        hess(:,iclose,:) = hess_near + hess_far;
    end
end

quasi_phase = exp(1i*kappa(:)*nx(:).'*d);


if nargout == 1
    val = quasi_phase.*val;

    if ising == 0
        for ii = -nsub:nsub
        isub = (abs(nx(:)-ii) > max(ls)) | ifar;

        if any(isub)
        vali = chnk.lap2d.green([0;0],[rx(isub).'+ (nx(isub).'-ii)*d;ry(isub).']);
        vali = reshape(vali,1,[],1);
        val(:,isub,:) = val(:,isub,:) - vali.*alpha.^(ii);
        end
        end
    end

    val = reshape(val,nkappa*ntarg,nsrc);
elseif nargout == 2
    val = quasi_phase.*val;
    grad = quasi_phase.*grad;
    
    if ising == 0
        for ii = -nsub:nsub        
        isub = (abs(nx(:)-ii) > max(ls)) | ifar;
    
        if any(isub)
        [vali, gradi] = chnk.lap2d.green([0;0],[rx(isub).' + (nx(isub).'-ii)*d;ry(isub).']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        val(:,isub,:) = val(:,isub,:) - vali.*alpha.^(ii);
        grad(:,isub,:) = grad(:,isub,:) - gradi.*alpha.^(ii);
        end
        end
    end

    val = reshape(val,nkappa*ntarg,nsrc);
    grad = reshape(grad,nkappa*ntarg,nsrc,2);
elseif nargout == 3
    val = quasi_phase.*val;
    grad = quasi_phase.*grad;
    hess = quasi_phase.*hess;
    
    if ising == 0
        for ii = -nsub:nsub        
        isub = (abs(nx(:)-ii) > max(ls)) | ifar;

        if any(isub)
        [vali, gradi, hessi] = chnk.lap2d.green([0;0],[rx(isub).' + (nx(isub).'-ii)*d;ry(isub).']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        hessi = reshape(hessi,1,[],3);

        val(:,isub,:) = val(:,isub,:) - vali.*alpha.^(ii);
        grad(:,isub,:) = grad(:,isub,:) - gradi.*alpha.^(ii);
        hess(:,isub,:) = hess(:,isub,:) - hessi.*alpha.^(ii);
        end
        end
    end

    val = reshape(val,nkappa*ntarg,nsrc);
    grad = reshape(grad,nkappa*ntarg,nsrc,2);
    hess = reshape(hess,nkappa*ntarg,nsrc,3);
end
end

% function [val, grad, hess] = green(src,targ,kappa,d,pxys,cs,l,ising)
% 
% pgreen = @(s,t) chnk.lap2d.green(s,t);
% sgreen = pgreen;
% 
% 
% [~,nsrc] = size(src);
% [~,ntarg] = size(targ);
% 
% xs = repmat(src(1,:),ntarg,1);
% ys = repmat(src(2,:),ntarg,1);
% 
% xt = repmat(targ(1,:).',1,nsrc);
% yt = repmat(targ(2,:).',1,nsrc);
% 
% rx = xt-xs;
% ry = yt-ys;
% 
% 
% rx = rx(:);
% ry = ry(:);
% 
% 
% nx = fix(rx/d);
% rx = rx - nx*d;
% 
% 
% rx2 = rx.*rx;
% ry2 = ry.*ry;
% 
% r2 = rx2+ry2;
% 
% r = sqrt(r2);
% 
% npt = size(r,1);
% 
% ythresh = 2*d/2;
% iclose = abs(ry) < ythresh;
% ifar = ~iclose;
% 
% rxfar = rx(ifar);
% ryfar = ry(ifar);
% 
% rxclose = rx(iclose);
% ryclose = ry(iclose);
% rclose = r(iclose);
% 
% nxclose = nx(iclose);
% nptclose = size(rxclose, 1);
% 
% nkappa = length(kappa);
% 
% val = zeros(nkappa,npt,1);
% if nargout > 1
% grad = zeros(nkappa,npt,2);
% end
% if nargout > 2
% hess = zeros(nkappa,npt,3);
% end
% 
% 
% if ~isempty(ryfar)
% 
%     if nargout > 2
%     [val(:,ifar,:), grad(:,ifar,:), hess(:,ifar,:)] =  quasi_dual_sum(rxfar,ryfar,0,kappa,d);
%     elseif nargout > 1
%     [val(:,ifar,:), grad(:,ifar,:)] =  quasi_dual_sum(rxfar,ryfar,0,kappa,d);
%     else
%     val(:,ifar,:) =  quasi_dual_sum(rxfar,ryfar,0,kappa,d);
%     end
% end
% 
% 
% alpha = (exp(1i*kappa(:)*d));
% 
% val_near= zeros(nkappa,nptclose);
% grad_near = zeros(nkappa,nptclose,2);
% hess_near = zeros(nkappa,nptclose,3);
% ls = -l:l;
% if ~isempty(rxclose)
%     for i = ls
%         if ising == 1
%             iuse = true(nptclose,1);
%         else
%             iuse = nxclose ~= -i;
%         end
% 
%         rxi = rxclose - i*d;
%         if nargout>2
%         [vali,gradi,hessi] = sgreen([0;0],[rxi.';ryclose.']);
%         vali = reshape(vali,1,[],1);
%         gradi = reshape(gradi,1,[],2);
%         hessi = reshape(hessi,1,[],3);
%         val_near(:,iuse) = val_near(:,iuse) + vali(:,iuse).*alpha.^i;
%         grad_near(:,iuse,:) = grad_near(:,iuse,:) + gradi(:,iuse,:).*alpha.^i;
%         hess_near(:,iuse,:) = hess_near(:,iuse,:) + hessi(:,iuse,:).*alpha.^i;
%         elseif nargout > 1
%         [vali,gradi] = sgreen([0;0],[rxi.';ryclose.']);
%         vali = reshape(vali,1,[],1);
%         gradi = reshape(gradi,1,[],2);
%         val_near(:,iuse) = val_near(:,iuse) + vali(:,iuse).*alpha.^i;
%         grad_near(:,iuse,:) = grad_near(:,iuse,:) + gradi(:,iuse,:).*alpha.^i;
%         else
%         vali = sgreen([0;0],[rxi.';ryclose.']);
%         vali = reshape(vali,1,[],1);
%         val_near(:,iuse) = val_near(:,iuse) + vali(:,iuse).*alpha.^i;
%         end
%     end
% 
% 
%     if nargout > 2
%         [val_far, grad_far, hess_far] = pgreen(pxys,[rxclose.';ryclose.']);
%         val_far = reshape(val_far * cs, nptclose, nkappa, 1);
%         grad_far = reshape(pagemtimes(grad_far, cs), nptclose, nkappa, 2);
%         hess_far = reshape(pagemtimes(hess_far, cs), nptclose, nkappa, 3);
% 
%         val(:,iclose,:) = val_near + permute(val_far, [2,1,3]); 
%         grad(:,iclose,:) = grad_near + permute(grad_far, [2,1,3]); 
%         hess(:,iclose,:) = hess_near + permute(hess_far, [2,1,3]);
%     elseif nargout >1
%         [val_far, grad_far] = pgreen(pxys,[rxclose.';ryclose.']);
%         val_far = reshape(val_far * cs, nptclose, nkappa, 1);
%         grad_far = reshape(pagemtimes(grad_far, cs), nptclose, nkappa, 2);
% 
%         val(:,iclose,:) = val_near + permute(val_far, [2,1,3]); 
%         grad(:,iclose,:) = grad_near + permute(grad_far, [2,1,3]);
%     else
%         val_far = pgreen(pxys,[rxclose.';ryclose.']);
%         val_far = reshape(val_far * cs, nptclose, nkappa, 1);
% 
%         val(:,iclose,:) = val_near + permute(val_far, [2,1,3]); 
%     end
% end
% 
% quasi_phase = exp(1i*kappa(:)*nx(:).'*d);
% 
% 
% if nargout == 1
%     val = quasi_phase.*val;
% 
%     if ising == 0
%         isub = (abs(nx(:)) > max(ls)) | ifar;
% 
%         if any(isub)
%         vali = sgreen([0;0],[rx(isub).'+ nx(isub).'*d;ry(isub).']);
%         vali = reshape(vali,1,[],1);
%         val(:,isub,:) = val(:,isub,:) - vali;
%         end
%     end
% 
%     val = reshape(val,nkappa*ntarg,nsrc);
% elseif nargout == 2
%     val = quasi_phase.*val;
%     grad = quasi_phase.*grad;
% 
%     if ising == 0
%         isub = (abs(nx(:)) > max(ls)) | ifar;
% 
%         if any(isub)
%         [vali, gradi] = sgreen([0;0],[rx(isub).' + nx(isub).'*d;ry(isub).']);
%         vali = reshape(vali,1,[],1);
%         gradi = reshape(gradi,1,[],2);
%         val(:,isub,:) = val(:,isub,:) - vali;
%         grad(:,isub,:) = grad(:,isub,:) - gradi;
%         end
%     end
% 
%     val = reshape(val,nkappa*ntarg,nsrc);
%     grad = reshape(grad,nkappa*ntarg,nsrc,2);
% elseif nargout == 3
%     val = quasi_phase.*val;
%     grad = quasi_phase.*grad;
%     hess = quasi_phase.*hess;
% 
%     if ising == 0
%         isub = (abs(nx(:)) > max(ls)) | ifar;
% 
%         if any(isub)
%         [vali, gradi, hessi] = sgreen([0;0],[rx(isub).' + nx(isub).'*d;ry(isub).']);
%         vali = reshape(vali,1,[],1);
%         gradi = reshape(gradi,1,[],2);
%         hessi = reshape(hessi,1,[],3);
% 
%         val(:,isub,:) = val(:,isub,:) - vali;
%         grad(:,isub,:) = grad(:,isub,:) - gradi;
%         hess(:,isub,:) = hess(:,isub,:) - hessi;
%         end
%     end
% 
%     val = reshape(val,nkappa*ntarg,nsrc);
%     grad = reshape(grad,nkappa*ntarg,nsrc,2);
%     hess = reshape(hess,nkappa*ntarg,nsrc,3);
% end
% 
% end
