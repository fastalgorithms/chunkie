function [val,grad,hess,third,fourth,fifth,sixth] = green(src,targ,zk,kappa,d,sn,l,ising,nsub)
%CHNK.FLEX2DQUAS.GREEN evaluate the quasiperiodic flexural wave Green's
% function for the given sources and targets.
%
% The quasiperiodic flexural Green's function is the difference of the
% quasiperiodic Helmholtz Green's functions at wavenumbers zk and i*zk:
%
%   G_flex(x,y) = G_H(x,y; zk) - G_H(x,y; i*zk)
%
% where G_H is the quasiperiodic Helmholtz Green's function satisfying
%   G_H(x + d e_1, y; k) = G_H(x, y; k) exp(i kappa d).
%
% G_flex is the fundamental solution to (Delta + zk^2)(Delta - zk^2) u = delta
% subject to the quasiperiodic condition u(x + d e_1) = u(x) exp(i kappa d).
%
% Syntax: [val,...] = chnk.flex2dquas.green(src,targ,zk,kappa,d,sn,l,ising)
%         [val,...] = chnk.flex2dquas.green(src,targ,zk,kappa,d,sn,l,ising,nsub)
%
% Input:
%   src   - (2,:) array of source positions
%   targ  - (2,:) array of target positions
%   zk    - complex number, flexural wavenumber
%   kappa - (nkappa,1) array of quasiperiodic phase parameters
%   d     - period (scalar)
%   sn    - (nkappa, N+1, 2) precomputed lattice sum coefficients;
%               sn(:,:,1) for the Helmholtz part at wavenumber zk,
%               sn(:,:,2) for the modified part at wavenumber i*zk
%               (see chnk.flex2dquas.latticecoefs)
%   l     - number of periodic copies included explicitly on each side
%   ising - if 1, include the free-space (singular) part; if 0, include
%               only the periodic images
%   nsub  - (optional, default 0) number of additional source copies to
%               subtract near the source
%
% Output:
%   val   - (nkappa*ntarg, nsrc) Green's function values
%   grad  - (nkappa*ntarg, nsrc, 2) gradient [d/dx, d/dy]
%   hess  - (nkappa*ntarg, nsrc, 3) Hessian [xx, xy, yy]
%   third - (nkappa*ntarg, nsrc, 4) third derivatives [xxx, xxy, xyy, yyy]
%   fourth- (nkappa*ntarg, nsrc, 5) fourth derivatives [xxxx, ..., yyyy]
%   fifth - (nkappa*ntarg, nsrc, 6) fifth derivatives
%   sixth - (nkappa*ntarg, nsrc, 7) sixth derivatives
%
% see also CHNK.FLEX2DQUAS.KERN, CHNK.FLEX2DQUAS.LATTICECOEFS,
%          CHNK.FLEX2DQUAS.HKDIFFGREEN1
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

ythresh = 2*d/2;
% ythresh = 2*d;
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

if nargin < 9
    nsub = 0;
elseif nsub > l
    error('trying to subtract off too many copies')
end

val = zeros(nkappa,npt,1);
if nargout > 1
grad = zeros(nkappa,npt,2);
end
if nargout > 2
hess = zeros(nkappa,npt,3);
end
if nargout > 3
third = zeros(nkappa,npt,4);
end
if nargout > 4
fourth = zeros(nkappa,npt,5);
end
if nargout > 5
fifth = zeros(nkappa,npt,6);
end
if nargout > 6
sixth = zeros(nkappa,npt,7);
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
betah = sqrt(1i*(xi_m-zk)).*sqrt(-1i*(xi_m+zk));
betak = sqrt(1i*(xi_m-1i*zk)).*sqrt(-1i*(xi_m+1i*zk));

fhath = exp(-betah.*sqrt(ryfar.^2) + 1i*xi_m.*rxfar)./(2*betah);
fhatk = exp(-betak.*sqrt(ryfar.^2) + 1i*xi_m.*rxfar)./(2*betak);

xifac = 1i*xi_m;
bfach = -betah.*(sqrt(ryfar.^2)./ryfar);
bfack = -betak.*(sqrt(ryfar.^2)./ryfar);

val(:,ifar,:) = sum(fhath - fhatk,3)/(d);
if nargout > 1
grad(:,ifar,1) = sum(xifac.*(fhath-fhatk),3)/d;
grad(:,ifar,2) = sum(bfach.*fhath-bfack.*fhatk,3)/d;
end

if nargout >2
hess(:,ifar,1) = sum(xifac.^2.*(fhath-fhatk),3)/d;
hess(:,ifar,2) = sum(xifac.*(bfach.*fhath-bfack.*fhatk),3)/d;
hess(:,ifar,3) = sum(bfach.^2.*fhath-bfack.^2.*fhatk,3)/d;
end

if nargout > 3
third(:,ifar,1) = sum(xifac.^3.*(fhath-fhatk),3)/d;
third(:,ifar,2) = sum(xifac.^2.*(bfach.*fhath-bfack.*fhatk),3)/d;
third(:,ifar,3) = sum(xifac.*(bfach.^2.*fhath-bfack.^2.*fhatk),3)/d;
third(:,ifar,4) = sum(bfach.^3.*fhath-bfack.^3.*fhatk,3)/d;
end

if nargout > 4
fourth(:,ifar,1) = sum(xifac.^4.*(fhath-fhatk),3)/d;
fourth(:,ifar,2) = sum(xifac.^3.*(bfach.*fhath-bfack.*fhatk),3)/d;
fourth(:,ifar,3) = sum(xifac.^2.*(bfach.^2.*fhath-bfack.^2.*fhatk),3)/d;
fourth(:,ifar,4) = sum(xifac.*(bfach.^3.*fhath-bfack.^3.*fhatk),3)/d;
fourth(:,ifar,5) = sum(bfach.^4.*fhath-bfack.^4.*fhatk,3)/d;
end

if nargout > 5
fifth(:,ifar,1) = sum(xifac.^5.*(fhath-fhatk),3)/d;
fifth(:,ifar,2) = sum(xifac.^4.*(bfach.*fhath-bfack.*fhatk),3)/d;
fifth(:,ifar,3) = sum(xifac.^3.*(bfach.^2.*fhath-bfack.^2.*fhatk),3)/d;
fifth(:,ifar,4) = sum(xifac.^2.*(bfach.^3.*fhath-bfack.^3.*fhatk),3)/d;
fifth(:,ifar,5) = sum(xifac.*(bfach.^4.*fhath-bfack.^4.*fhatk),3)/d;
fifth(:,ifar,6) = sum(bfach.^5.*fhath-bfack.^5.*fhatk,3)/d;
end

if nargout > 6
sixth(:,ifar,1) = sum(xifac.^6.*(fhath-fhatk),3)/d;
sixth(:,ifar,2) = sum(xifac.^5.*(bfach.*fhath-bfack.*fhatk),3)/d;
sixth(:,ifar,3) = sum(xifac.^4.*(bfach.^2.*fhath-bfack.^2.*fhatk),3)/d;
sixth(:,ifar,4) = sum(xifac.^3.*(bfach.^3.*fhath-bfack.^3.*fhatk),3)/d;
sixth(:,ifar,5) = sum(xifac.^2.*(bfach.^4.*fhath-bfack.^4.*fhatk),3)/d;
sixth(:,ifar,6) = sum(xifac.*(bfach.^5.*fhath-bfack.^5.*fhatk),3)/d;
sixth(:,ifar,7) = sum(bfach.^6.*fhath-bfack.^6.*fhatk,3)/d;
end

end

alpha = (exp(1i*kappa(:)*d));

val_near= zeros(nkappa,nptclose);
grad_near = zeros(nkappa,nptclose,2);
hess_near = zeros(nkappa,nptclose,3);
third_near = zeros(nkappa,nptclose,4);
fourth_near = zeros(nkappa,nptclose,5);
fifth_near = zeros(nkappa,nptclose,6);
sixth_near = zeros(nkappa,nptclose,7);
ls = -l:l;
if ~isempty(rxclose)
    for i = ls
        if ising == 1
            iuse = true(nptclose,1);
        else
            iuse = ~ismember(nxclose, -i-nsub:-i+nsub); 
        end

        rxi = rxclose - i*d;
        if nargout>6
        [vali,gradi,hessi,thirdi,fourthi,fifthi,sixthi] = chnk.flex2dquas.hkdiffgreen1(zk,[0;0],[rxi.';ryclose.']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        hessi = reshape(hessi,1,[],3);
        thirdi = reshape(thirdi,1,[],4);
        fourthi = reshape(fourthi,1,[],5);
        fifthi = reshape(fifthi,1,[],6);
        sixthi = reshape(sixthi,1,[],7);
        val_near(:,iuse) = val_near(:,iuse) + vali(:,iuse).*alpha.^i;
        grad_near(:,iuse,:) = grad_near(:,iuse,:) + gradi(:,iuse,:).*alpha.^i;
        hess_near(:,iuse,:) = hess_near(:,iuse,:) + hessi(:,iuse,:).*alpha.^i;   
        third_near(:,iuse,:) = third_near(:,iuse,:) + thirdi(:,iuse,:).*alpha.^i;   
        fourth_near(:,iuse,:) = fourth_near(:,iuse,:) + fourthi(:,iuse,:).*alpha.^i;   
        fifth_near(:,iuse,:) = fifth_near(:,iuse,:) + fifthi(:,iuse,:).*alpha.^i;   
        sixth_near(:,iuse,:) = sixth_near(:,iuse,:) + sixthi(:,iuse,:).*alpha.^i;   
        elseif nargout>5
        [vali,gradi,hessi,thirdi,fourthi,fifthi] = chnk.flex2d.hkdiffgreen(zk,[0;0],[rxi.';ryclose.']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        hessi = reshape(hessi,1,[],3);
        thirdi = reshape(thirdi,1,[],4);
        fourthi = reshape(fourthi,1,[],5);
        fifthi = reshape(fifthi,1,[],6);
        val_near(:,iuse) = val_near(:,iuse) + vali(:,iuse).*alpha.^i;
        grad_near(:,iuse,:) = grad_near(:,iuse,:) + gradi(:,iuse,:).*alpha.^i;
        hess_near(:,iuse,:) = hess_near(:,iuse,:) + hessi(:,iuse,:).*alpha.^i;   
        third_near(:,iuse,:) = third_near(:,iuse,:) + thirdi(:,iuse,:).*alpha.^i;   
        fourth_near(:,iuse,:) = fourth_near(:,iuse,:) + fourthi(:,iuse,:).*alpha.^i;   
        fifth_near(:,iuse,:) = fifth_near(:,iuse,:) + fifthi(:,iuse,:).*alpha.^i;   
        elseif nargout>4
        [vali,gradi,hessi,thirdi,fourthi] = chnk.flex2d.hkdiffgreen(zk,[0;0],[rxi.';ryclose.']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        hessi = reshape(hessi,1,[],3);
        thirdi = reshape(thirdi,1,[],4);
        fourthi = reshape(fourthi,1,[],5);
        val_near(:,iuse) = val_near(:,iuse) + vali(:,iuse).*alpha.^i;
        grad_near(:,iuse,:) = grad_near(:,iuse,:) + gradi(:,iuse,:).*alpha.^i;
        hess_near(:,iuse,:) = hess_near(:,iuse,:) + hessi(:,iuse,:).*alpha.^i;   
        third_near(:,iuse,:) = third_near(:,iuse,:) + thirdi(:,iuse,:).*alpha.^i;   
        fourth_near(:,iuse,:) = fourth_near(:,iuse,:) + fourthi(:,iuse,:).*alpha.^i;   
        elseif nargout>3
        [vali,gradi,hessi,thirdi] = chnk.flex2d.hkdiffgreen(zk,[0;0],[rxi.';ryclose.']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        hessi = reshape(hessi,1,[],3);
        thirdi = reshape(thirdi,1,[],4);
        val_near(:,iuse) = val_near(:,iuse) + vali(:,iuse).*alpha.^i;
        grad_near(:,iuse,:) = grad_near(:,iuse,:) + gradi(:,iuse,:).*alpha.^i;
        hess_near(:,iuse,:) = hess_near(:,iuse,:) + hessi(:,iuse,:).*alpha.^i;   
        third_near(:,iuse,:) = third_near(:,iuse,:) + thirdi(:,iuse,:).*alpha.^i;        
        elseif nargout>2
        [vali,gradi,hessi] = chnk.flex2d.hkdiffgreen(zk,[0;0],[rxi.';ryclose.']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        hessi = reshape(hessi,1,[],3);
        val_near(:,iuse) = val_near(:,iuse) + vali(:,iuse).*alpha.^i;
        grad_near(:,iuse,:) = grad_near(:,iuse,:) + gradi(:,iuse,:).*alpha.^i;
        hess_near(:,iuse,:) = hess_near(:,iuse,:) + hessi(:,iuse,:).*alpha.^i;
        elseif nargout > 1
        [vali,gradi] = chnk.flex2d.hkdiffgreen(zk,[0;0],[rxi.';ryclose.']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        val_near(:,iuse) = val_near(:,iuse) + vali(:,iuse).*alpha.^i;
        grad_near(:,iuse,:) = grad_near(:,iuse,:) + gradi(:,iuse,:).*alpha.^i;
        else
        vali = chnk.flex2d.hkdiffgreen(zk,[0;0],[rxi.';ryclose.']);
        vali = reshape(vali,1,[],1);
        val_near(:,iuse) = val_near(:,iuse) + vali(:,iuse).*alpha.^i;
        end
    end

    N = size(sn,2)-1;
    ns = (0:N);
    ns_use = (0:N+7);
    Js = zeros(length(rclose),N+8);
    Is = zeros(length(rclose),N+8);

    if length(rclose) < N+8
        for i = 1:length(rclose)
        Js(i,:) = besselj(ns_use,zk*rclose(i));
        Is(i,:) = besselj(ns_use,1i*zk*rclose(i));
        end
    else
        for i = 1:length(ns_use)
        Js(:,i) = besselj(ns_use(i),zk*rclose);
        Is(:,i) = besselj(ns_use(i),1i*zk*rclose);
        end
    end
    % t1 = tic;
    eip = (rxclose+1i*ryclose)./rclose;
    eipn = reshape(eip.^ns,1,[], N+1);
    cs = (eipn+1./eipn)/2;
    % cs(isnan(cs)) = 0;
    
    Js = reshape(Js,1,[],N+8);
    Is = reshape(Is,1,[],N+8);
    snj = reshape(sn(:,:,1),nkappa, 1, N+1);
    sni = reshape(sn(:,:,2),nkappa, 1, N+1);
    
    tmpj = reshape(Js(:,:,2:end-7).*cs(:,:,2:end),[],N);
    tmpi = reshape(-Is(:,:,2:end-7).*cs(:,:,2:end),[],N);    
    val_far = 0.25*1i*(Js(:,:,1).*snj(:,:,1)-Is(:,:,1).*sni(:,:,1)) + 0.5*1i*snj(:,2:end)*tmpj.' + 0.5*1i*sni(:,2:end)*tmpi.';
    
    ndiag = sum(rclose < 1e-14);
    if ndiag > 0
    val_far(:,rclose < 1e-14) = repmat(0.25*1i*(snj(:,:,1)-sni(:,:,1)),1,ndiag); % diagonal replacement
    end
    val(:,iclose) = val_near+val_far;
    
    if nargout >1
        DJs = cat(3,-Js(:,:,2),.5*(Js(:,:,1:end-8)-Js(:,:,3:end-6)))*zk;
        DIs = cat(3,-Is(:,:,2),.5*(Is(:,:,1:end-8)-Is(:,:,3:end-6)))*1i*zk;
        ss = (eipn-1./eipn)/2i;
        % ss(isnan(ss)) = 0;
        
        tmpj = reshape(DJs(:,:,2:end).*cs(:,:,2:end),[],N);
        tmpi = reshape(-DIs(:,:,2:end).*cs(:,:,2:end),[],N);

        grad_far_p = 0.25*1i*(DJs(:,:,1).*snj(:,:,1)-DIs(:,:,1).*sni(:,:,1)) + 0.5*1i*snj(:,2:end)*tmpj.' + 0.5*1i*sni(:,2:end)*tmpi.';

        tmpj = reshape((Js(:,:,2:end-7)).*ss(:,:,2:end),[],N)./rclose;
        tmpi = reshape((-Is(:,:,2:end-7)).*ss(:,:,2:end),[],N)./rclose;

        grad_far_t = (0.5*1i*((-reshape((1:N),1,[]).*sni(:,2:end))*tmpi.'))+(0.5*1i*((-reshape((1:N),1,[]).*snj(:,2:end))*tmpj.'));
        
        grad_far = cat(3,cs(:,:,2).*grad_far_p - ss(:,:,2).*grad_far_t, ss(:,:,2).*grad_far_p + cs(:,:,2).*grad_far_t);

        if ndiag > 0
        grad_far(:,rclose < 1e-14,1) = repmat(0.25*1i*(zk*snj(:,:,2)-1i*zk*sni(:,:,2)),1,ndiag);
        grad_far(:,rclose < 1e-14,2) = 0;
        end

        grad(:,iclose,:) = grad_near + grad_far; 
    end
    if nargout > 2
        DDJs = cat(3,.5*(Js(:,:,3)-Js(:,:,1)),.25*(Js(:,:,4)-3*Js(:,:,2)),.25*(Js(:,:,1:end-9)-2*Js(:,:,3:end-7)+Js(:,:,5:end-5)))*zk^2;
        DDIs = cat(3,.5*(Is(:,:,3)-Is(:,:,1)),.25*(Is(:,:,4)-3*Is(:,:,2)),.25*(Is(:,:,1:end-9)-2*Is(:,:,3:end-7)+Is(:,:,5:end-5)))*(1i*zk)^2;
        rclose = rclose.';
        rxclose = rxclose.';
        ryclose = ryclose.';
        ns = reshape(ns,1,1,[]);

        tmp_nj = rclose.^(-4).*(-ns.*ryclose.*(Js(:,:,1:end-7)).*(ns.*ryclose.*cs+2*rxclose.*ss)+ ...
            rclose.*ryclose.*(ryclose.*cs + 2*ns.*rxclose.*ss).*(DJs)+ ...
            rclose.^2.*rxclose.^2.*cs.*(DDJs));
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = rclose.^(-4).*(-ns.*ryclose.*(-Is(:,:,1:end-7)).*(ns.*ryclose.*cs+2*rxclose.*ss)+ ...
            rclose.*ryclose.*(ryclose.*cs + 2*ns.*rxclose.*ss).*(-DIs)+ ...
            rclose.^2.*rxclose.^2.*cs.*(-DDIs));
        tmp_ni = reshape(tmp_ni,[],N+1);
        hess_far_xx = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';

        tmp_nj = ns.*(Js(:,:,1:end-7)).*(ns.*rxclose.*ryclose.*cs+(rxclose.^2-ryclose.^2).*ss).*rclose.^(-4) ...
            +rclose.^(-3).*(-(rxclose.*ryclose.*cs+ns.*(rxclose.^2-ryclose.^2).*ss).*(DJs)+...
            rclose.*rxclose.*ryclose.*cs.*(DDJs));
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = ns.*(-Is(:,:,1:end-7)).*(ns.*rxclose.*ryclose.*cs+(rxclose.^2-ryclose.^2).*ss).*rclose.^(-4) ...
            +rclose.^(-3).*(-(rxclose.*ryclose.*cs+ns.*(rxclose.^2-ryclose.^2).*ss).*(-DIs)+...
            rclose.*rxclose.*ryclose.*cs.*(-DDIs));
        tmp_ni = reshape(tmp_ni,[],N+1);
        hess_far_xy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';

        tmp_nj = rclose.^(-4).*(-ns.*rxclose.*(Js(:,:,1:end-7)).*(ns.*rxclose.*cs-2*ryclose.*ss)+ ...
            rclose.*rxclose.*(rxclose.*cs - 2*ns.*ryclose.*ss).*(DJs)+ ...
            rclose.^2.*ryclose.^2.*cs.*(DDJs));
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = rclose.^(-4).*(-ns.*rxclose.*(-Is(:,:,1:end-7)).*(ns.*rxclose.*cs-2*ryclose.*ss)+ ...
            rclose.*rxclose.*(rxclose.*cs - 2*ns.*ryclose.*ss).*(-DIs)+ ...
            rclose.^2.*ryclose.^2.*cs.*(-DDIs));
        tmp_ni = reshape(tmp_ni,[],N+1);
        hess_far_yy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';

        hess_far = cat(3,hess_far_xx, hess_far_xy, hess_far_yy);

        if ndiag > 0     
        hess_far(:,rclose < 1e-14,1) = repmat(-0.25*1i*zk^2/2*(snj(:,:,1)+sni(:,:,1))+.5*1i*zk^2/4*(snj(:,:,3)+sni(:,:,3)),1,ndiag);
        hess_far(:,rclose < 1e-14,2) = 0;
        hess_far(:,rclose < 1e-14,3) = repmat(-0.25*1i*zk^2/2*(snj(:,:,1)+sni(:,:,1))-.5*1i*zk^2/4*(snj(:,:,3)+sni(:,:,3)),1,ndiag);
        end

        hess(:,iclose,:) = hess_near + hess_far;
    end
    if nargout > 3
        D3Js = cat(3,1/4*(3*Js(:,:,2) - Js(:,:,4)), ...
            1/8*(4*Js(:,:,3)-Js(:,:,5)-3*Js(:,:,1)), ...
            1/8*(-4*Js(:,:,2)+3*Js(:,:,4)-Js(:,:,6)), ...
            1/8*(Js(:,:,1:end-10)-3*Js(:,:,3:end-8)+3*Js(:,:,5:end-6)-Js(:,:,7:end-4))) *zk^3;
        D3Is = cat(3,1/4*(3*Is(:,:,2) - Is(:,:,4)), ...
            1/8*(4*Is(:,:,3)-Is(:,:,5)-3*Is(:,:,1)), ...
            1/8*(-4*Is(:,:,2)+3*Is(:,:,4)-Is(:,:,6)), ...
            1/8*(Is(:,:,1:end-10)-3*Is(:,:,3:end-8)+3*Is(:,:,5:end-6)-Is(:,:,7:end-4))) *(1i*zk)^3;

        tmp_nj = rclose.^(-6).*(ns.*ryclose.*(Js(:,:,1:end-7)).*(6*ns.*rxclose.*ryclose.*cs+(6*rxclose.^2 - (2+ns.^2).*ryclose.^2).*ss)+ ...
            rclose.*(-3.*ryclose.*((1+ns.^2).*rxclose.*ryclose.*cs + ns.*(2*rxclose.^2-ryclose.^2).*ss)).*(DJs)+ ...
            rclose.^2.*3.*rxclose.*ryclose.*(ryclose.*cs+ns.*rxclose.*ss).*(DDJs) + ...
            rclose.^3.*rxclose.^3.*cs.*(D3Js));
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = rclose.^(-6).*(ns.*ryclose.*(-Is(:,:,1:end-7)).*(6*ns.*rxclose.*ryclose.*cs+(6*rxclose.^2 - (2+ns.^2).*ryclose.^2).*ss)+ ...
            rclose.*(-3.*ryclose.*((1+ns.^2).*rxclose.*ryclose.*cs + ns.*(2*rxclose.^2-ryclose.^2).*ss)).*(-DIs)+ ...
            rclose.^2.*3.*rxclose.*ryclose.*(ryclose.*cs+ns.*rxclose.*ss).*(-DDIs) + ...
            rclose.^3.*rxclose.^3.*cs.*(-D3Is));
        tmp_ni = reshape(tmp_ni,[],N+1);
        third_far_xxx = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';


        tmp_nj = rclose.^(-6).*(ns.*(Js(:,:,1:end-7)).*(2*ns.*ryclose.*(-2*rxclose.^2+ryclose.^2).*cs + rxclose.*(-2*rxclose.^2+(6+ns.^2).*ryclose.^2).*ss) + ...
            rclose.*(-(1+ns.^2).*ryclose.*(-2*rxclose.^2+ryclose.^2).*cs + ns.*rxclose.*(2*rxclose.^2-7*ryclose.^2).*ss).*(DJs) + ...
            -rclose.^2.*(-ryclose.*(-2*rxclose.^2 + ryclose.^2).*cs + ns.*rxclose.*(rxclose.^2 - 2*ryclose.^2).*ss).*(DDJs) + ...
            rclose.^3.*(rxclose.^2.*ryclose).*cs.*(D3Js)); 
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = rclose.^(-6).*(ns.*(-Is(:,:,1:end-7)).*(2*ns.*ryclose.*(-2*rxclose.^2+ryclose.^2).*cs + rxclose.*(-2*rxclose.^2+(6+ns.^2).*ryclose.^2).*ss) + ...
            rclose.*(-(1+ns.^2).*ryclose.*(-2*rxclose.^2+ryclose.^2).*cs + ns.*rxclose.*(2*rxclose.^2-7*ryclose.^2).*ss).*(-DIs) + ...
            -rclose.^2.*(-ryclose.*(-2*rxclose.^2 + ryclose.^2).*cs + ns.*rxclose.*(rxclose.^2 - 2*ryclose.^2).*ss).*(-DDIs) + ...
            rclose.^3.*(rxclose.^2.*ryclose).*cs.*(-D3Is)); 
        tmp_ni = reshape(tmp_ni,[],N+1);
        third_far_xxy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';

        tmp_nj = rclose.^(-6).*(ns.*(Js(:,:,1:end-7)).*(2*ns.*rxclose.*(-2*ryclose.^2+rxclose.^2).*cs + ryclose.*(2*ryclose.^2-(6+ns.^2).*rxclose.^2).*ss) + ...
            rclose.*(-(1+ns.^2).*rxclose.*(-2*ryclose.^2+rxclose.^2).*cs - ns.*ryclose.*(2*ryclose.^2-7*rxclose.^2).*ss).*(DJs) + ...
            rclose.^2.*(rxclose.*(-2*ryclose.^2 + rxclose.^2).*cs + ns.*ryclose.*(ryclose.^2 - 2*rxclose.^2).*ss).*(DDJs) + ...
            rclose.^3.*(ryclose.^2.*rxclose).*cs.*(D3Js)); 
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = rclose.^(-6).*(ns.*(-Is(:,:,1:end-7)).*(2*ns.*rxclose.*(-2*ryclose.^2+rxclose.^2).*cs + ryclose.*(2*ryclose.^2-(6+ns.^2).*rxclose.^2).*ss) + ...
            rclose.*(-(1+ns.^2).*rxclose.*(-2*ryclose.^2+rxclose.^2).*cs - ns.*ryclose.*(2*ryclose.^2-7*rxclose.^2).*ss).*(-DIs) + ...
            rclose.^2.*(rxclose.*(-2*ryclose.^2 + rxclose.^2).*cs + ns.*ryclose.*(ryclose.^2 - 2*rxclose.^2).*ss).*(-DDIs) + ...
            rclose.^3.*(ryclose.^2.*rxclose).*cs.*(-D3Is)); 
        tmp_ni = reshape(tmp_ni,[],N+1);
        third_far_xyy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';

        tmp_nj = rclose.^(-6).*(ns.*rxclose.*(Js(:,:,1:end-7)).*(6*ns.*ryclose.*rxclose.*cs+(-6*ryclose.^2 + (2+ns.^2).*rxclose.^2).*ss)+ ...
            rclose.*(-3.*rxclose.*((1+ns.^2).*ryclose.*rxclose.*cs + ns.*(-2*ryclose.^2+rxclose.^2).*ss)).*(DJs)+ ...
            rclose.^2.*3.*ryclose.*rxclose.*(rxclose.*cs-ns.*ryclose.*ss).*(DDJs) + ...
            rclose.^3.*ryclose.^3.*cs.*(D3Js));
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = rclose.^(-6).*(ns.*rxclose.*(-Is(:,:,1:end-7)).*(6*ns.*ryclose.*rxclose.*cs+(-6*ryclose.^2 + (2+ns.^2).*rxclose.^2).*ss)+ ...
            rclose.*(-3.*rxclose.*((1+ns.^2).*ryclose.*rxclose.*cs + ns.*(-2*ryclose.^2+rxclose.^2).*ss)).*(-DIs)+ ...
            rclose.^2.*3.*ryclose.*rxclose.*(rxclose.*cs-ns.*ryclose.*ss).*(-DDIs) + ...
            rclose.^3.*ryclose.^3.*cs.*(-D3Is));
        tmp_ni = reshape(tmp_ni,[],N+1);
        third_far_yyy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';

        third_far = cat(3,third_far_xxx,third_far_xxy,third_far_xyy,third_far_yyy);

        if ndiag > 0      
        third_far(:,rclose < 1e-14,1) = repmat(-.5*1i*(6*(zk)^3.*snj(:,:,2)-6*(1i*zk)^3*sni(:,:,2))/16 + .5*1i*(6*(zk)^3.*snj(:,:,4)-6*(1i*zk)^3*sni(:,:,4))/48,1,ndiag) ;
        third_far(:,rclose < 1e-14,2) = 0;
        third_far(:,rclose < 1e-14,3) = repmat(-.5*1i*(2*(zk)^3.*snj(:,:,2)-2*(1i*zk)^3*sni(:,:,2))/16 - .5*1i*(6*(zk)^3.*snj(:,:,4)-6*(1i*zk)^3*sni(:,:,4))/48,1,ndiag) ;
        third_far(:,rclose < 1e-14,4) = 0;
        end

        third(:,iclose,:) = third_near + third_far;
    end
    if nargout > 4
        D4Js = cat(3,1/8*(3*Js(:,:,1)-4*Js(:,:,3)+Js(:,:,5)), ...
            1/16*(10*Js(:,:,2)-5*Js(:,:,4)+Js(:,:,6)), ...
            1/16*(-4*Js(:,:,1)+7*Js(:,:,3)-4*Js(:,:,5)+Js(:,:,7)   ), ...
            1/16*(-5*Js(:,:,2)+6*Js(:,:,4)-4*Js(:,:,6)+Js(:,:,8)   ), ...
            1/16*(Js(:,:,1:end-11)-4*Js(:,:,3:end-9)+6*Js(:,:,5:end-7)-4*Js(:,:,7:end-5)+Js(:,:,9:end-3))    ) *zk^4;
        D4Is = cat(3,1/8*(3*Is(:,:,1)-4*Is(:,:,3)+Is(:,:,5)), ...
            1/16*(10*Is(:,:,2)-5*Is(:,:,4)+Is(:,:,6)), ...
            1/16*(-4*Is(:,:,1)+7*Is(:,:,3)-4*Is(:,:,5)+Is(:,:,7)   ), ...
            1/16*(-5*Is(:,:,2)+6*Is(:,:,4)-4*Is(:,:,6)+Is(:,:,8)   ), ...
            1/16*(Is(:,:,1:end-11)-4*Is(:,:,3:end-9)+6*Is(:,:,5:end-7)-4*Is(:,:,7:end-5)+Is(:,:,9:end-3))   ) *zk^4;

        tmp_nj = rclose.^(-8).*(ns.*ryclose.*(Js(:,:,1:end-7)).*(ns.*ryclose.*(-36*rxclose.^2+(8+ns.^2).*ryclose.^2).*cs + 12*rxclose.*(-2*rxclose.^2+(2+ns.^2).*ryclose.^2).*ss) + ...
            rclose.*(ryclose.*(-3*(1+2*ns.^2).*ryclose.*(-4*rxclose.^2+ryclose.^2).*cs - 4*ns.*rxclose.*(-6*rxclose.^2 + (8+ns.^2).*ryclose.^2).*ss).*(DJs)) + ...
            rclose.^2.*(3*ryclose.*((-2*(2+ns.^2).*rxclose.^2.*ryclose + ryclose.^3).*cs + 4*ns.*rxclose.*(-rxclose.^2 + ryclose.^2).*ss)).*(DDJs) + ...
            rclose.^3.*rxclose.^2.*(2*ryclose.*(3*ryclose.*cs + 2*ns.*rxclose.*ss)).*(D3Js) + ...
            rclose.^4.*rxclose.^4.*cs.*(D4Js));
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = rclose.^(-8).*(ns.*ryclose.*(-Is(:,:,1:end-7)).*(ns.*ryclose.*(-36*rxclose.^2+(8+ns.^2).*ryclose.^2).*cs + 12*rxclose.*(-2*rxclose.^2+(2+ns.^2).*ryclose.^2).*ss) + ...
            rclose.*(ryclose.*(-3*(1+2*ns.^2).*ryclose.*(-4*rxclose.^2+ryclose.^2).*cs - 4*ns.*rxclose.*(-6*rxclose.^2 + (8+ns.^2).*ryclose.^2).*ss).*(-DIs)) + ...
            rclose.^2.*(3*ryclose.*((-2*(2+ns.^2).*rxclose.^2.*ryclose + ryclose.^3).*cs + 4*ns.*rxclose.*(-rxclose.^2 + ryclose.^2).*ss)).*(-DDIs) + ...
            rclose.^3.*rxclose.^2.*(2*ryclose.*(3*ryclose.*cs + 2*ns.*rxclose.*ss)).*(-D3Is) + ...
            rclose.^4.*rxclose.^4.*cs.*(-D4Is));
        tmp_ni = reshape(tmp_ni,[],N+1);
        fourth_far_xxxx = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';

        tmp_nj = rclose.^(-8).*(ns.*(Js(:,:,1:end-7)).*(-ns.*rxclose.*ryclose.*(-18*rxclose.^2 + (26+ns.^2).*ryclose.^2).*cs + 3*(2*rxclose.^4 - 3*(4+ns.^2).*rxclose.^2.*ryclose.^2+(2+ns.^2).*ryclose.^4).*ss) ...
            - rclose.*(3*(1+2*ns.^2).*rxclose.*ryclose.*(2*rxclose.^2-3*ryclose.^2).*cs + ns.*(6*rxclose.^4-3*(14+ns.^2).*rxclose.^2.*ryclose.^2 + (8+ns.^2).*ryclose.^4).*ss).*(DJs) ...
            + 3*rclose.^2.*(rxclose.*ryclose.*((2+ns.^2).*rxclose.^2 - (3+ns.^2).*ryclose.^2).*cs + ns.*(rxclose.^4 - 6*rxclose.^2.*ryclose.^2 + ryclose.^4).*ss).*(DDJs) - ...
            rxclose.*rclose.^3.*((3*ryclose.*(rxclose.^2-ryclose.^2).*cs + ns.*rxclose.*(rxclose.^2 - 3*ryclose.^2).*ss).*(D3Js) - rxclose.^2.*ryclose.*rclose.*cs.*(D4Js)));
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = rclose.^(-8).*(ns.*(-Is(:,:,1:end-7)).*(-ns.*rxclose.*ryclose.*(-18*rxclose.^2 + (26+ns.^2).*ryclose.^2).*cs + 3*(2*rxclose.^4 - 3*(4+ns.^2).*rxclose.^2.*ryclose.^2+(2+ns.^2).*ryclose.^4).*ss) ...
            - rclose.*(3*(1+2*ns.^2).*rxclose.*ryclose.*(2*rxclose.^2-3*ryclose.^2).*cs + ns.*(6*rxclose.^4-3*(14+ns.^2).*rxclose.^2.*ryclose.^2 + (8+ns.^2).*ryclose.^4).*ss).*(-DIs) ...
            + 3*rclose.^2.*(rxclose.*ryclose.*((2+ns.^2).*rxclose.^2 - (3+ns.^2).*ryclose.^2).*cs + ns.*(rxclose.^4 - 6*rxclose.^2.*ryclose.^2 + ryclose.^4).*ss).*(-DDIs) - ...
            rxclose.*rclose.^3.*((3*ryclose.*(rxclose.^2-ryclose.^2).*cs + ns.*rxclose.*(rxclose.^2 - 3*ryclose.^2).*ss).*(-D3Is) - rxclose.^2.*ryclose.*rclose.*cs.*(-D4Is)));
        tmp_ni = reshape(tmp_ni,[],N+1);
        fourth_far_xxxy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';

        tmp_nj = rclose.^(-8).*(ns.*(Js(:,:,1:end-7)).*(ns.*(-6*rxclose.^4 + (32+ns.^2).*rxclose.^2.*ryclose.^2 - 6*ryclose.^4).*cs + 6*(4+ns.^2).*rxclose.*ryclose.*(rxclose.^2 - ryclose.^2).*ss) + ...
            rclose.*((1+2*ns.^2).*(2*rxclose.^4 - 11*rxclose.^2.*ryclose.^2 + 2*ryclose.^4).*cs - 2*ns.*(14+ns.^2).*rxclose.*ryclose.*(rxclose.^2 - ryclose.^2).*ss).*(DJs) - ...
            rclose.^2.*(((2+ns.^2).*rxclose.^4 - (11+4*ns.^2).*rxclose.^2.*ryclose.^2 + (2+ns.^2).*ryclose.^4).*cs + 12*ns.*rxclose.*ryclose.*(-rxclose.^2 + ryclose.^2).*ss).*(DDJs) + ...
            rclose.^3.*((rxclose.^4 - 4*rxclose.^2.*ryclose.^2 + ryclose.^4).*cs + 2*ns.*rxclose.*ryclose.*(-rxclose.^2 + ryclose.^2).*ss).*(D3Js) + ...
            rclose.^4.*(rxclose.^2.*ryclose.^2.*cs).*(D4Js));
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = rclose.^(-8).*(ns.*(-Is(:,:,1:end-7)).*(ns.*(-6*rxclose.^4 + (32+ns.^2).*rxclose.^2.*ryclose.^2 - 6*ryclose.^4).*cs + 6*(4+ns.^2).*rxclose.*ryclose.*(rxclose.^2 - ryclose.^2).*ss) + ...
            rclose.*((1+2*ns.^2).*(2*rxclose.^4 - 11*rxclose.^2.*ryclose.^2 + 2*ryclose.^4).*cs - 2*ns.*(14+ns.^2).*rxclose.*ryclose.*(rxclose.^2 - ryclose.^2).*ss).*(-DIs) - ...
            rclose.^2.*(((2+ns.^2).*rxclose.^4 - (11+4*ns.^2).*rxclose.^2.*ryclose.^2 + (2+ns.^2).*ryclose.^4).*cs + 12*ns.*rxclose.*ryclose.*(-rxclose.^2 + ryclose.^2).*ss).*(-DDIs) + ...
            rclose.^3.*((rxclose.^4 - 4*rxclose.^2.*ryclose.^2 + ryclose.^4).*cs + 2*ns.*rxclose.*ryclose.*(-rxclose.^2 + ryclose.^2).*ss).*( - D3Is) + ...
            rclose.^4.*(rxclose.^2.*ryclose.^2.*cs).*(-D4Is));
        tmp_ni = reshape(tmp_ni,[],N+1);
        fourth_far_xxyy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';

        tmp_nj = rclose.^(-8).*(-ns.*(Js(:,:,1:end-7)).*(ns.*rxclose.*ryclose.*(-18*ryclose.^2 + (26+ns.^2).*rxclose.^2).*cs + 3*(2*ryclose.^4 - 3*(4+ns.^2).*rxclose.^2.*ryclose.^2+(2+ns.^2).*rxclose.^4).*ss) ...
            + rclose.*(3*(1+2*ns.^2).*rxclose.*ryclose.*(-2*ryclose.^2+3*rxclose.^2).*cs + ns.*(6*ryclose.^4-3*(14+ns.^2).*rxclose.^2.*ryclose.^2 + (8+ns.^2).*rxclose.^4).*ss).*(DJs) ...
            - 3*rclose.^2.*(rxclose.*ryclose.*(-(2+ns.^2).*ryclose.^2 + (3+ns.^2).*rxclose.^2).*cs + ns.*(rxclose.^4 - 6*rxclose.^2.*ryclose.^2 + ryclose.^4).*ss).*(DDJs) + ...
            ryclose.*rclose.^3.*((3*rxclose.*(rxclose.^2-ryclose.^2).*cs + ns.*ryclose.*(ryclose.^2 - 3*rxclose.^2).*ss).*(D3Js) + ryclose.^2.*rxclose.*rclose.*cs.*(D4Js)));
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = rclose.^(-8).*(-ns.*(-Is(:,:,1:end-7)).*(ns.*rxclose.*ryclose.*(-18*ryclose.^2 + (26+ns.^2).*rxclose.^2).*cs + 3*(2*ryclose.^4 - 3*(4+ns.^2).*rxclose.^2.*ryclose.^2+(2+ns.^2).*rxclose.^4).*ss) ...
            + rclose.*(3*(1+2*ns.^2).*rxclose.*ryclose.*(-2*ryclose.^2+3*rxclose.^2).*cs + ns.*(6*ryclose.^4-3*(14+ns.^2).*rxclose.^2.*ryclose.^2 + (8+ns.^2).*rxclose.^4).*ss).*(-DIs) ...
            - 3*rclose.^2.*(rxclose.*ryclose.*(-(2+ns.^2).*ryclose.^2 + (3+ns.^2).*rxclose.^2).*cs + ns.*(rxclose.^4 - 6*rxclose.^2.*ryclose.^2 + ryclose.^4).*ss).*(-DDIs) + ...
            ryclose.*rclose.^3.*((3*rxclose.*(rxclose.^2-ryclose.^2).*cs + ns.*ryclose.*(ryclose.^2 - 3*rxclose.^2).*ss).*(-D3Is) + ryclose.^2.*rxclose.*rclose.*cs.*(-D4Is)));
        tmp_ni = reshape(tmp_ni,[],N+1);
        fourth_far_xyyy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';

        tmp_nj = rclose.^(-8).*(ns.*rxclose.*(Js(:,:,1:end-7)).*(ns.*rxclose.*(-36*ryclose.^2+(8+ns.^2).*rxclose.^2).*cs - 12*ryclose.*(-2*ryclose.^2+(2+ns.^2).*rxclose.^2).*ss) + ...
            rclose.*(-rxclose.*(3*(1+2*ns.^2).*rxclose.*(-4*ryclose.^2+rxclose.^2).*cs - 4*ns.*ryclose.*(-6*ryclose.^2 + (8+ns.^2).*rxclose.^2).*ss).*(DJs)) + ...
            rclose.^2.*(3*rxclose.*((-2*(2+ns.^2).*ryclose.^2.*rxclose + rxclose.^3).*cs + 4*ns.*ryclose.*(-rxclose.^2 + ryclose.^2).*ss)).*(DDJs) + ...
            rclose.^3.*ryclose.^2.*(2*rxclose.*(3*rxclose.*cs - 2*ns.*ryclose.*ss)).*(D3Js) + ...
            rclose.^4.*ryclose.^4.*cs.*(D4Js));
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = rclose.^(-8).*(ns.*rxclose.*(-Is(:,:,1:end-7)).*(ns.*rxclose.*(-36*ryclose.^2+(8+ns.^2).*rxclose.^2).*cs - 12*ryclose.*(-2*ryclose.^2+(2+ns.^2).*rxclose.^2).*ss) + ...
            rclose.*(-rxclose.*(3*(1+2*ns.^2).*rxclose.*(-4*ryclose.^2+rxclose.^2).*cs - 4*ns.*ryclose.*(-6*ryclose.^2 + (8+ns.^2).*rxclose.^2).*ss).*(-DIs)) + ...
            rclose.^2.*(3*rxclose.*((-2*(2+ns.^2).*ryclose.^2.*rxclose + rxclose.^3).*cs + 4*ns.*ryclose.*(-rxclose.^2 + ryclose.^2).*ss)).*(-DDIs) + ...
            rclose.^3.*ryclose.^2.*(2*rxclose.*(3*rxclose.*cs - 2*ns.*ryclose.*ss)).*(-D3Is) + ...
            rclose.^4.*ryclose.^4.*cs.*(-D4Is));
        tmp_ni = reshape(tmp_ni,[],N+1);
        fourth_far_yyyy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';        

        fourth_far = cat(3,fourth_far_xxxx,fourth_far_xxxy,fourth_far_xxyy,fourth_far_xyyy,fourth_far_yyyy);

        if ndiag > 0        
        fourth_far(:,rclose < 1e-14,1) = repmat(0.25*1i*(24*(zk)^4.*snj(:,:,1)-24*(1i*zk)^4.*sni(:,:,1))/(2^6) - 0.5*1i*(24*(zk)^4.*snj(:,:,3)-24*(1i*zk)^4.*sni(:,:,3))/(2^4*6) + 0.5*1i*(24*(zk)^4.*snj(:,:,5)-24*(1i*zk)^4.*sni(:,:,5))/(2^4*24),1,ndiag) ;
        fourth_far(:,rclose < 1e-14,2) = 0;
        fourth_far(:,rclose < 1e-14,3) = repmat(0.25*1i*(8*(zk)^4.*snj(:,:,1)-8*(1i*zk)^4.*sni(:,:,1))/(2^6) + 0*0.5*1i*(24*(zk)^4.*snj(:,:,3)-24*(1i*zk)^4.*sni(:,:,3))/(2^4*6) - 0.5*1i*(24*(zk)^4.*snj(:,:,5)-24*(1i*zk)^4.*sni(:,:,5))/(2^4*24),1,ndiag) ;
        fourth_far(:,rclose < 1e-14,4) = 0;
        fourth_far(:,rclose < 1e-14,5) = repmat(0.25*1i*(24*(zk)^4.*snj(:,:,1)-24*(1i*zk)^4.*sni(:,:,1))/(2^6) + 0.5*1i*(24*(zk)^4.*snj(:,:,3)-24*(1i*zk)^4.*sni(:,:,3))/(2^4*6) + 0.5*1i*(24*(zk)^4.*snj(:,:,5)-24*(1i*zk)^4.*sni(:,:,5))/(2^4*24),1,ndiag) ;
        end

        fourth(:,iclose,:) = fourth_near + fourth_far;
    end
    if nargout > 5
        D5Js = cat(3,1/16*(-10*Js(:,:,2)+5*Js(:,:,4)-Js(:,:,6)), ...
            1/32*(10*Js(:,:,1)-15*Js(:,:,3)+6*Js(:,:,5)-Js(:,:,7)), ...
            1/32*(15*Js(:,:,2)-11*Js(:,:,4)+5*Js(:,:,6)-Js(:,:,8)   ), ...
            1/32*(-5*Js(:,:,1)+11*Js(:,:,3)-10*Js(:,:,5)+5*Js(:,:,7)-Js(:,:,9)   ), ...
            1/32*(-6*Js(:,:,2)+10*Js(:,:,4)-10*Js(:,:,6)+5*Js(:,:,8)-Js(:,:,10)),  ...
            1/32*(Js(:,:,1:end-12)-5*Js(:,:,3:end-10) +10*Js(:,:,5:end-8)-10*Js(:,:,7:end-6)+5*Js(:,:,9:end-4)-Js(:,:,11:end-2))) *zk^5;
        D5Is = cat(3,1/16*(-10*Is(:,:,2)+5*Is(:,:,4)-Is(:,:,6)), ...
            1/32*(10*Is(:,:,1)-15*Is(:,:,3)+6*Is(:,:,5)-Is(:,:,7)), ...
            1/32*(15*Is(:,:,2)-11*Is(:,:,4)+5*Is(:,:,6)-Is(:,:,8)   ), ...
            1/32*(-5*Is(:,:,1)+11*Is(:,:,3)-10*Is(:,:,5)+5*Is(:,:,7)-Is(:,:,9)   ), ...
            1/32*(-6*Is(:,:,2)+10*Is(:,:,4)-10*Is(:,:,6)+5*Is(:,:,8)-Is(:,:,10)),  ...
            1/32*(Is(:,:,1:end-12)-5*Is(:,:,3:end-10) +10*Is(:,:,5:end-8)-10*Is(:,:,7:end-6)+5*Is(:,:,9:end-4)-Is(:,:,11:end-2))) *(1i*zk^5);

        x = rxclose; x2 = x.*x; x3 = x2.*x; x4 = x2.*x2;
        y = ryclose; y2 = y.*y; y3 = y2.*y; y4 = y2.*y2;
        n = ns; n2 = n.*n; n4 = n2.*n2;
        r2 = rclose.^2;
        r9 = rclose.^9;
        r10 = r2.^5;

        fj = Js(:,:,1:end-7);
        fjp = DJs;
        fjpp = DDJs;
        fjp3 = D3Js;
        fjp4 = D4Js;
        fjp5 = D5Js;

        fi = -Is(:,:,1:end-7);
        fip = -DIs;
        fipp = -DDIs;
        fip3 = -D3Is;
        fip4 = -D4Is;
        fip5 = -D5Is;

        tmp_nj = 1./r10.*(n.*y.*fj.*(-20.*n.*x.*y.*(-12*x2+(8+n2).*y2).*cs + (120*x4 - 120*(2+n2).*x2.*y2 + (24 + 20*n2+n4).*y4).*ss ) + ...
            rclose.*(-5*y.*(x.*y.*(12*(1+3*n2).*x2-(9+26*n2+n4).*y2).*cs +n.*(24*x4-12*(5+n2).*x2.*y2+(7+2*n2).*y4).*ss).*fjp + ...
            5*y.*rclose.*(3*x.*y.*(4*(1+n2).*x2-(3+2*n2).*y2).*cs+n.*(12*x4-2*(14+n2).*x2.*y2+3*y4).*ss).*fjpp+ ... 
            x.*r2.*(5*y.*(y.*(-2*(3+n2).*x2+3*y2).*cs+2*n.*x.*(-2*x2+3 *y2).* ss).* fjp3+5 *x2.* y .* rclose .* (2 * y .* cs+n .* x .* ss).*fjp4 +x4.*r2.* cs .* fjp5)));
        tmp_nj = reshape(tmp_nj,[],N+1);

        tmp_ni = (n.*y.*fi.*(-20.*n.*x.*y.*(-12*x2+(8+n2).*y2).*cs + (120*x4 - 120*(2+n2).*x2.*y2 + (24 + 20*n2+n4).*y4).*ss ))./r10 + ...
            1./r9.*(-5*y.*(x.*y.*(12*(1+3*n2).*x2-(9+26*n2+n4).*y2).*cs +n.*(24*x4-12*(5+n2).*x2.*y2+(7+2*n2).*y4).*ss).*fip + ...
            5*y.*rclose.*(3*x.*y.*(4*(1+n2).*x2-(3+2*n2).*y2).*cs+n.*(12*x4-2*(14+n2).*x2.*y2+3*y4).*ss).*fipp+ ... 
            x.*r2.*(5*y.*(y.*(-2*(3+n2).*x2+3*y2).*cs+2*n.*x.*(-2*x2+3 *y2).* ss).* fip3+5 *x2.* y .* rclose .* (2 * y .* cs+n .* x .* ss).*fip4 +x4.* (x2+y2).* cs .* fip5));
        tmp_ni = reshape(tmp_ni,[],N+1);

        fifth_far_xxxxx = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';

        tmp_nj = -n.*fj.*(4.*n.*y.*(24*x4-4*(17+n2).*x2.*y2+(8+n2).*y4).*cs + x.*(24*x4-24*(10+3*n2).*x2.*y2+(120+68*n2+n4).*y4).*ss)./r10 ...
                + 1./r9.*( (y.*(24*(1+3*n2).*x4-4*(18+53*n2+n4).*x2.*y2+(9+26*n2+n4).*y4).*cs ...
                + n.*x.*(24*x4-12*(23+3*n2).*x2.*y2+(155+34*n2).*y4).*ss).*fjp ...
                - rclose.*(3*y.*(8*(1+n2).*x4-4*(6+5*n2).*x2.*y2+(3+2*n2).*y4).*cs ...
                + n.*x.*(12*x4-6*(22+n2).*x2.*y2+(71+4*n2).*y4).*ss).*fjpp ...
                + r2.*((y.*(4*(3+n2).*x4-6*(5+n2).*x2.*y2+3*y4).*cs ...
                + 2*n.*x.*(2*x4-17*x2.*y2+6*y4).*ss).*fjp3 ...
                - x2.*rclose.*((4*x2.*y-6*y3).*cs + n.*x.*(x2-4*y2).*ss).*fjp4 + x4.*y.*r2.*cs.*fjp5 ) );
        tmp_nj = reshape(tmp_nj,[],N+1);

        tmp_ni = -n.*fi.*(4.*n.*y.*(24*x4-4*(17+n2).*x2.*y2+(8+n2).*y4).*cs + x.*(24*x4-24*(10+3*n2).*x2.*y2+(120+68*n2+n4).*y4).*ss)./r10 ...
                + 1./r9.*( (y.*(24*(1+3*n2).*x4-4*(18+53*n2+n4).*x2.*y2+(9+26*n2+n4).*y4).*cs ...
                + n.*x.*(24*x4-12*(23+3*n2).*x2.*y2+(155+34*n2).*y4).*ss).*fip ...
                - rclose.*(3*y.*(8*(1+n2).*x4-4*(6+5*n2).*x2.*y2+(3+2*n2).*y4).*cs ...
                + n.*x.*(12*x4-6*(22+n2).*x2.*y2+(71+4*n2).*y4).*ss).*fipp ...
                + r2.*((y.*(4*(3+n2).*x4-6*(5+n2).*x2.*y2+3*y4).*cs ...
                + 2*n.*x.*(2*x4-17*x2.*y2+6*y4).*ss).*fip3 ...
                - x2.*rclose.*((4*x2.*y-6*y3).*cs + n.*x.*(x2-4*y2).*ss).*fip4 + x4.*y.*r2.*cs.*fip5 ) );
        tmp_ni = reshape(tmp_ni,[],N+1);

        fifth_far_xxxxy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';

        tmp_nj = (n.*fj.*(-4.*n.*x.*(-6*x4+3*(20+n2).*x2.*y2-2*(17+n2).*y4).*cs ...
                   + y.*(-12*(10+3*n2).*x4+(240+92*n2+n4).*x2.*y2-12*(2+n2).*y4).*ss))./r10 ...
                   + 1./r9.*((-x.*(6*(1+3*n2).*x4-3*(21+62*n2+n4).*x2.*y2+2*(18+53*n2+n4).*y4).*cs ...
                   + n.*y.*(6*(23+3*n2).*x4-(287+46*n2).*x2.*y2+6*(5+n2).*y4).*ss).*fjp ...
                   + rclose.*(3.*x.*(2*(1+n2).*x4-3*(7+6*n2).*x2.*y2+2*(6+5*n2).*y4).*cs ...
                   - n.*y.*(3*(22+n2).*x4-3*(45+2*n2).*x2.*y2+(14+n2).*y4).*ss).*fjpp ...
                   + r2.*( -(x.*((3+n2).*x4-3*(9+2*n2).*x2.*y2+3*(5+n2).*y4).*cs ...
                   + n.*y.*(-17*x4+30*x2.*y2-3*y4).*ss).*fjp3 + x.*rclose.*((x4-6*x2.*y2+3*y4).*cs ...
                   + n.*x.*y.*(-2*x2+3*y2).*ss).*fjp4 + x3.*y2.*r2.*cs.*fjp5));
        tmp_nj = reshape(tmp_nj,[],N+1);

        tmp_ni = (n.*fi.*(-4.*n.*x.*(-6*x4+3*(20+n2).*x2.*y2-2*(17+n2).*y4).*cs ...
                   + y.*(-12*(10+3*n2).*x4+(240+92*n2+n4).*x2.*y2-12*(2+n2).*y4).*ss))./r10 ...
                   + 1./r9.*((-x.*(6*(1+3*n2).*x4-3*(21+62*n2+n4).*x2.*y2+2*(18+53*n2+n4).*y4).*cs ...
                   + n.*y.*(6*(23+3*n2).*x4-(287+46*n2).*x2.*y2+6*(5+n2).*y4).*ss).*fip ...
                   + rclose.*(3.*x.*(2*(1+n2).*x4-3*(7+6*n2).*x2.*y2+2*(6+5*n2).*y4).*cs ...
                   - n.*y.*(3*(22+n2).*x4-3*(45+2*n2).*x2.*y2+(14+n2).*y4).*ss).*fipp ...
                   + r2.*( -(x.*((3+n2).*x4-3*(9+2*n2).*x2.*y2+3*(5+n2).*y4).*cs ...
                   + n.*y.*(-17*x4+30*x2.*y2-3*y4).*ss).*fip3 + x.*rclose.*((x4-6*x2.*y2+3*y4).*cs ...
                   + n.*x.*y.*(-2*x2+3*y2).*ss).*fip4 + x3.*y2.*r2.*cs.*fip5));
        tmp_ni = reshape(tmp_ni,[],N+1);

        fifth_far_xxxyy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';   

        tmp_nj = n.*fj.*(4.*n.*y.*(2*(17+n2).*x4-3*(20+n2).*x2.*y2+6*y4).*cs ...
          + x.*(12*(2+n2).*x4-(240+92*n2+n4).*x2.*y2+12*(10+3*n2).*y4).*ss)./r10 ...
          - 1./r9.*( (y.*(2*(18+53*n2+n4).*x4-3*(21+62*n2+n4).*x2.*y2+6*(1+3*n2).*y4).*cs ...
          + n.*x.*(6*(5+n2).*x4-(287+46*n2).*x2.*y2+6*(23+3*n2).*y4).*ss).*fjp ...
          - rclose.*(3.*y.*(2*(6+5*n2).*x4-3*(7+6*n2).*x2.*y2+2*(1+n2).*y4).*cs ...
            + n.*x.*((14+n2).*x4-3*(45+2*n2).*x2.*y2+3*(22+n2).*y4).*ss).*fjpp ...
          + r2.*( (y.*(3*(5+n2).*x4-3*(9+2*n2).*x2.*y2+(3+n2).*y4).*cs ...
          + n.*x.*(3*x4-30*x2.*y2+17*y4).*ss).*fjp3 ...
          - y.*( rclose.*((3*x4-6*x2.*y2+y4).*cs + n.*x.*y.*(-3*x2+2*y2).*ss).*fjp4 ...
          + x2.*y2.*r2.*cs.*fjp5 ) ) );
        tmp_nj = reshape(tmp_nj,[],N+1);

        tmp_ni = n.*fi.*(4.*n.*y.*(2*(17+n2).*x4-3*(20+n2).*x2.*y2+6*y4).*cs ...
          + x.*(12*(2+n2).*x4-(240+92*n2+n4).*x2.*y2+12*(10+3*n2).*y4).*ss)./r10 ...
          - 1./r9.*( (y.*(2*(18+53*n2+n4).*x4-3*(21+62*n2+n4).*x2.*y2+6*(1+3*n2).*y4).*cs ...
          + n.*x.*(6*(5+n2).*x4-(287+46*n2).*x2.*y2+6*(23+3*n2).*y4).*ss).*fip ...
          - rclose.*(3.*y.*(2*(6+5*n2).*x4-3*(7+6*n2).*x2.*y2+2*(1+n2).*y4).*cs ...
            + n.*x.*((14+n2).*x4-3*(45+2*n2).*x2.*y2+3*(22+n2).*y4).*ss).*fipp ...
          + r2.*( (y.*(3*(5+n2).*x4-3*(9+2*n2).*x2.*y2+(3+n2).*y4).*cs ...
          + n.*x.*(3*x4-30*x2.*y2+17*y4).*ss).*fip3 ...
          - y.*( rclose.*((3*x4-6*x2.*y2+y4).*cs + n.*x.*y.*(-3*x2+2*y2).*ss).*fip4 ...
          + x2.*y2.*r2.*cs.*fip5 ) ) );
        tmp_ni = reshape(tmp_ni,[],N+1);

        fifth_far_xxyyy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';   

        tmp_nj = n.*fj.*(-4.*n.*x.*((8+n2).*x4-4*(17+n2).*x2.*y2+24*y4).*cs ...
            + y.*((120+68*n2+n4).*x4-24*(10+3*n2).*x2.*y2+24*y4).*ss)./r10 ...
            + 1./r9.*( (x.*((9+26*n2+n4).*x4-4*(18+53*n2+n4).*x2.*y2+24*(1+3*n2).*y4).*cs ...
            + n.*y.*(-(155+34*n2).*x4+12*(23+3*n2).*x2.*y2-24*y4).*ss).*fjp ...
            - rclose.*(3.*x.*((3+2*n2).*x4-4*(6+5*n2).*x2.*y2+8*(1+n2).*y4).*cs ...
            + n.*y.*(-(71+4*n2).*x4+6*(22+n2).*x2.*y2-12*y4).*ss).*fjpp ...
            + r2.*( (x.*(3*x4-6*(5+n2).*x2.*y2+4*(3+n2).*y4).*cs ...
            - 2.*n.*y.*(6*x4-17*x2.*y2+2*y4).*ss).*fjp3 ...
            + y2.*rclose.*((6*x3-4*x.*y2).*cs + n.*y.*(-4*x2+y2).*ss).*fjp4 ...
            + x.*y4.*r2.*cs.*fjp5 ) );
        tmp_nj = reshape(tmp_nj,[],N+1);

        tmp_ni = n.*fi.*(-4.*n.*x.*((8+n2).*x4-4*(17+n2).*x2.*y2+24*y4).*cs ...
            + y.*((120+68*n2+n4).*x4-24*(10+3*n2).*x2.*y2+24*y4).*ss)./r10 ...
            + 1./r9.*( (x.*((9+26*n2+n4).*x4-4*(18+53*n2+n4).*x2.*y2+24*(1+3*n2).*y4).*cs ...
            + n.*y.*(-(155+34*n2).*x4+12*(23+3*n2).*x2.*y2-24*y4).*ss).*fip ...
            - rclose.*(3.*x.*((3+2*n2).*x4-4*(6+5*n2).*x2.*y2+8*(1+n2).*y4).*cs ...
            + n.*y.*(-(71+4*n2).*x4+6*(22+n2).*x2.*y2-12*y4).*ss).*fipp ...
            + r2.*( (x.*(3*x4-6*(5+n2).*x2.*y2+4*(3+n2).*y4).*cs ...
            - 2.*n.*y.*(6*x4-17*x2.*y2+2*y4).*ss).*fip3 ...
            + y2.*rclose.*((6*x3-4*x.*y2).*cs + n.*y.*(-4*x2+y2).*ss).*fip4 ...
            + x.*y4.*r2.*cs.*fip5 ) );
        tmp_ni = reshape(tmp_ni,[],N+1);

        fifth_far_xyyyy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';   

        tmp_nj = -n.*x.*fj.*(20.*n.*x.*y.*((8+n2).*x2-12.*y2).*cs + ((24+20*n2+n4).*x4-120.*(2+n2).*x2.*y2+120.*y4).*ss)./r10 ...
            + 1./r9.*(  5.*x.*(x.*y.*((9+26*n2+n4).*x2-12.*(1+3*n2).*y2).*cs ...
            + n.*((7+2*n2).*x4-12.*(5+n2).*x2.*y2+24.*y4).*ss).*fjp ...
            - 5.*x.*rclose.*(3.*x.*y.*((3+2*n2).*x2-4.*(1+n2).*y2).*cs ...
            + n.*(3.*x4-2.*(14+n2).*x2.*y2+12.*y4).*ss).*fjpp ...
            + y.*r2.*(  5.*x.*(x.*(3.*x2-2.*(3+n2).*y2).*cs + 2.*n.*y.*(-3.*x2+2.*y2).*ss).*fjp3 ...
            + 5.*x.*y2.*rclose.*(2.*x.*cs - n.*y.*ss).*fjp4 + y4.*r2.*cs.*fjp5));
        tmp_nj = reshape(tmp_nj,[],N+1);
        
        tmp_ni = -n.*x.*fi.*(20.*n.*x.*y.*((8+n2).*x2-12.*y2).*cs + ((24+20*n2+n4).*x4-120.*(2+n2).*x2.*y2+120.*y4).*ss)./r10 ...
            + 1./r9.*(  5.*x.*(x.*y.*((9+26*n2+n4).*x2-12.*(1+3*n2).*y2).*cs ...
            + n.*((7+2*n2).*x4-12.*(5+n2).*x2.*y2+24.*y4).*ss).*fip ...
            - 5.*x.*rclose.*(3.*x.*y.*((3+2*n2).*x2-4.*(1+n2).*y2).*cs ...
            + n.*(3.*x4-2.*(14+n2).*x2.*y2+12.*y4).*ss).*fipp ...
            + y.*r2.*(  5.*x.*(x.*(3.*x2-2.*(3+n2).*y2).*cs + 2.*n.*y.*(-3.*x2+2.*y2).*ss).*fip3 ...
            + 5.*x.*y2.*rclose.*(2.*x.*cs - n.*y.*ss).*fip4 + y4.*r2.*cs.*fip5)); 
        tmp_ni = reshape(tmp_ni,[],N+1);

        fifth_far_yyyyy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';   

        fifth_far = cat(3,fifth_far_xxxxx,fifth_far_xxxxy,fifth_far_xxxyy,fifth_far_xxyyy,fifth_far_xyyyy,fifth_far_yyyyy);

        if ndiag > 0
        fifth_far(:,rclose < 1e-14,1) = repmat(.5*1i*(120*(zk)^5.*snj(:,:,2)-120*(1i*zk)^5*sni(:,:,2))/(6*2^6) - .5*1i*(120*(zk)^5.*snj(:,:,4)-120*(1i*zk)^5*sni(:,:,4))/(24*2^5) + .5*1i*(120*(zk)^5*snj(:,:,6)-120*(1i*zk)^5*sni(:,:,6))/(120*2^5),1,ndiag);
        fifth_far(:,rclose < 1e-14,2) = 0;
        fifth_far(:,rclose < 1e-14,3) = repmat(.5*1i*(24*(zk)^5.*snj(:,:,2)-24*(1i*zk)^5*sni(:,:,2))/(6*2^6) + .5*1i*(24*(zk)^5.*snj(:,:,4)-24*(1i*zk)^5*sni(:,:,4))/(24*2^5) - .5*1i*(120*(zk)^5*snj(:,:,6)-120*(1i*zk)^5*sni(:,:,6))/(120*2^5),1,ndiag);
        fifth_far(:,rclose < 1e-14,4) = 0;
        fifth_far(:,rclose < 1e-14,5) = repmat(.5*1i*(24*(zk)^5.*snj(:,:,2)-24*(1i*zk)^5*sni(:,:,2))/(6*2^6) + .5*1i*(72*(zk)^5.*snj(:,:,4)-72*(1i*zk)^5*sni(:,:,4))/(24*2^5) + .5*1i*(120*(zk)^5*snj(:,:,6)-120*(1i*zk)^5*sni(:,:,6))/(120*2^5),1,ndiag);
        fifth_far(:,rclose < 1e-14,6) = 0;
        end

        fifth(:,iclose,:) = fifth_near + fifth_far;
    end
    if nargout > 6
        D6Js = cat(3,1/32*(-10*Js(:,:,1)+15*Js(:,:,3)-6*Js(:,:,5)+Js(:,:,7)), ...
            1/64*(-35*Js(:,:,2)+21*Js(:,:,4)-7*Js(:,:,6)+Js(:,:,8)), ...
            1/64*(15*Js(:,:,1)-26*Js(:,:,3)+16*Js(:,:,5)-6*Js(:,:,7)+Js(:,:,9)), ...
            1/64*(21*Js(:,:,2)-21*Js(:,:,4)+15*Js(:,:,6)-6*Js(:,:,8)+Js(:,:,10)), ...
            1/64*(-6*Js(:,:,1)+16*Js(:,:,3)-20*Js(:,:,5)+15*Js(:,:,7)-6*Js(:,:,9)+Js(:,:,11)),  ...
            1/64*(-7*Js(:,:,2)+15*Js(:,:,4)-20*Js(:,:,6)+15*Js(:,:,8)-6*Js(:,:,10)+Js(:,:,12) ), ...
            1/64*(Js(:,:,1:end-13)-6*Js(:,:,3:end-11)+15*Js(:,:,5:end-9)-20*Js(:,:,7:end-7)+15*Js(:,:,9:end-5)-6*Js(:,:,11:end-3)+Js(:,:,13:end-1))) *zk^6;
        D6Is = cat(3,1/32*(-10*Is(:,:,1)+15*Is(:,:,3)-6*Is(:,:,5)+Is(:,:,7)), ...
            1/64*(-35*Is(:,:,2)+21*Is(:,:,4)-7*Is(:,:,6)+Is(:,:,8)), ...
            1/64*(15*Is(:,:,1)-26*Is(:,:,3)+16*Is(:,:,5)-6*Is(:,:,7)+Is(:,:,9)), ...
            1/64*(21*Is(:,:,2)-21*Is(:,:,4)+15*Is(:,:,6)-6*Is(:,:,8)+Is(:,:,10)), ...
            1/64*(-6*Is(:,:,1)+16*Is(:,:,3)-20*Is(:,:,5)+15*Is(:,:,7)-6*Is(:,:,9)+Is(:,:,11)),  ...
            1/64*(-7*Is(:,:,2)+15*Is(:,:,4)-20*Is(:,:,6)+15*Is(:,:,8)-6*Is(:,:,10)+Is(:,:,12) ), ...
            1/64*(Is(:,:,1:end-13)-6*Is(:,:,3:end-11)+15*Is(:,:,5:end-9)-20*Is(:,:,7:end-7)+15*Is(:,:,9:end-5)-6*Is(:,:,11:end-3)+Is(:,:,13:end-1))) *(1i*zk)^6;

        x6 = x4.*x2;
        y6 = y4.*y2;
        n4 = n2.*n2;
        r12 = r2.^6;

        fjp6 = D6Js;
        fip6 = -D6Is;

        tmp_nj = (-n.*y.*fj.*( n.*y.*(1800.*x4-300.*(8+n2).*x2.*y2+(184+40.*n2+n4).*y4).*cs ...
         + 30.*x.*(24.*x4-40.*(2+n2).*x2.*y2+(24+20.*n2+n4).*y4).*ss) + 3.*y.*rclose.*( ...
           5.*y.*(24.*(1+4.*n2).*x4-4.*(9+34.*n2+2.*n4).*x2.*y2+(3+11.*n2+n4).*y4).*cs ...
         + 2.*n.*x.*(120.*x4-120.*(4+n2).*x2.*y2+(164+60.*n2+n4).*y4).*ss ).*fjp ...
         - 15.*y.*r2.*( y.*(12.*(2+3.*n2).*x4-(36+44.*n2+n4).*x2.*y2+3.*(1+n2).*y4).*cs ...
         + 4.*n.*x.*(6.*x4-3.*(8+n2).*x2.*y2+(8+n2).*y4).*ss ).*fjpp + r2.*rclose.*( ...
           5.*y.*( 3.*y.*(4.*(3+2.*n2).*x4-3.*(5+2.*n2).*x2.*y2+y4).*cs ...
         + 2.*n.*x.*(12.*x4-2.*(20+n2).*x2.*y2+9.*y4).*ss ).*fjp3 ...
         - 15.*x2.*y.*rclose.*(y.*((4+n2).*x2-3.*y2).*cs + 2.*n.*x.*(x2-2.*y2).*ss).*fjp4 ...
         + x4.*r2.*(3.*y.*(5.*y.*cs + 2.*n.*x.*ss).*fjp5 + x2.*rclose.*cs.*fjp6)))./r12 ; 
        tmp_nj = reshape(tmp_nj,[],N+1);

        tmp_ni = (-n.*y.*fi.*( n.*y.*(1800.*x4-300.*(8+n2).*x2.*y2+(184+40.*n2+n4).*y4).*cs ...
         + 30.*x.*(24.*x4-40.*(2+n2).*x2.*y2+(24+20.*n2+n4).*y4).*ss) + 3.*y.*rclose.*( ...
          5.*y.*(24.*(1+4.*n2).*x4-4.*(9+34.*n2+2.*n4).*x2.*y2+(3+11.*n2+n4).*y4).*cs ...
         + 2.*n.*x.*(120.*x4-120.*(4+n2).*x2.*y2+(164+60.*n2+n4).*y4).*ss ).*fip ...
         - 15.*y.*r2.*( y.*(12.*(2+3.*n2).*x4-(36+44.*n2+n4).*x2.*y2+3.*(1+n2).*y4).*cs ...
         + 4.*n.*x.*(6.*x4-3.*(8+n2).*x2.*y2+(8+n2).*y4).*ss ).*fipp + r2.*rclose.*( ...
           5.*y.*( 3.*y.*(4.*(3+2.*n2).*x4-3.*(5+2.*n2).*x2.*y2+y4).*cs ...
         + 2.*n.*x.*(12.*x4-2.*(20+n2).*x2.*y2+9.*y4).*ss ).*fip3 ...
         - 15.*x2.*y.*rclose.*(y.*((4+n2).*x2-3.*y2).*cs + 2.*n.*x.*(x2-2.*y2).*ss).*fip4 ...
         + x4.*r2.*(3.*y.*(5.*y.*cs + 2.*n.*x.*ss).*fip5 + x2.*rclose.*cs.*fip6)))./r12 ;         
        tmp_ni = reshape(tmp_ni,[],N+1);

        sixth_far_xxxxxx = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';  

        tmp_nj = (n.*fj.*(n.*x.*y.*(600.*x4-200.*(14+n2).*x2.*y2+(984+140.*n2+n4).*y4).*cs ...
             + 5.*(24.*x6-120.*(3+n2).*x4.*y2+5.*(72+44.*n2+n4).*x2.*y4-(24+20.*n2+n4).*y6).*ss ) ...
             - rclose.*((5.*x.*y.*(24.*(1+4.*n2).*x4-8.*(15+58.*n2+2.*n4).*x2.*y2+(45+169.*n2+11.*n4).*y4).*cs ...
             + n.*(120.*x6-120.*(17+3.*n2).*x4.*y2+5.*(452+132.*n2+n4).*x2.*y4-(164+60.*n2+n4).*y6).*ss).*fjp ...
             - 5.*rclose.*(x.*y.*(12.*(2+3.*n2).*x4-2.*(60+80.*n2+n4).*x2.*y2+(45+53.*n2+n4).*y4).*cs ...
             - 2.*n.*(-6.*x6+3.*(34+3.*n2).*x4.*y2-14.*(8+n2).*x2.*y4+(8+n2).*y6).*ss).*fjpp ...
             + r2.*( 5.*(2.*x.*y.*((6+4.*n2).*x4-(27+14.*n2).*x2.*y2+3.*(3+n2).*y4).*cs ...
             + n.*(4.*x6-2.*(30+n2).*x4.*y2+(55+2.*n2).*x2.*y4-3.*y6).*ss).*fjp3 ...
             - 5.*x.*rclose.*(y.*((4+n2).*x4-2.*(7+n2).*x2.*y2+3.*y4).*cs ...
             + n.*x.*(x4-11.*x2.*y2+6.*y4).*ss).*fjp4 + x3.*r2.*( ...
              (5.*y.*(x2-2.*y2).*cs + n.*x.*(x2-5.*y2).*ss).*fjp5 - x2.*y.*rclose.*cs.*fjp6 ))))./r12;
        tmp_nj = reshape(tmp_nj,[],N+1);

        tmp_ni = (n.*fi.*(n.*x.*y.*(600.*x4-200.*(14+n2).*x2.*y2+(984+140.*n2+n4).*y4).*cs ...
             + 5.*(24.*x6-120.*(3+n2).*x4.*y2+5.*(72+44.*n2+n4).*x2.*y4-(24+20.*n2+n4).*y6).*ss ) ...
             - rclose.*((5.*x.*y.*(24.*(1+4.*n2).*x4-8.*(15+58.*n2+2.*n4).*x2.*y2+(45+169.*n2+11.*n4).*y4).*cs ...
             + n.*(120.*x6-120.*(17+3.*n2).*x4.*y2+5.*(452+132.*n2+n4).*x2.*y4-(164+60.*n2+n4).*y6).*ss).*fip ...
             - 5.*rclose.*(x.*y.*(12.*(2+3.*n2).*x4-2.*(60+80.*n2+n4).*x2.*y2+(45+53.*n2+n4).*y4).*cs ...
             - 2.*n.*(-6.*x6+3.*(34+3.*n2).*x4.*y2-14.*(8+n2).*x2.*y4+(8+n2).*y6).*ss).*fipp ...
             + r2.*( 5.*(2.*x.*y.*((6+4.*n2).*x4-(27+14.*n2).*x2.*y2+3.*(3+n2).*y4).*cs ...
             + n.*(4.*x6-2.*(30+n2).*x4.*y2+(55+2.*n2).*x2.*y4-3.*y6).*ss).*fip3 ...
             - 5.*x.*rclose.*(y.*((4+n2).*x4-2.*(7+n2).*x2.*y2+3.*y4).*cs ...
             + n.*x.*(x4-11.*x2.*y2+6.*y4).*ss).*fip4 + x3.*r2.*( ...
              (5.*y.*(x2-2.*y2).*cs + n.*x.*(x2-5.*y2).*ss).*fip5 - x2.*y.*rclose.*cs.*fip6 ))))./r12;        
        tmp_ni = reshape(tmp_ni,[],N+1);

        sixth_far_xxxxxy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';  

        tmp_nj = (-n.*fj.*(n.*(120.*x6-120.*(16+n2).*x4.*y2+(2184+200.*n2+n4).*x2.*y4-20.*(8+n2).*y6).*cs ...
  - 10.*x.*y.*(24.*(3+n2).*x4-2.*(120+56.*n2+n4).*x2.*y2+(72+44.*n2+n4).*y4).*ss) ...
  + rclose.*((24.*(1+4.*n2).*x6-24.*(17+66.*n2+2.*n4).*x4.*y2+(477+1829.*n2+79.*n4).*x2.*y4 ...
       - 4.*(9+34.*n2+2.*n4).*y6).*cs+ 2.*n.*x.*y.*(-24.*(17+3.*n2).*x4+2.*(716+168.*n2+n4).*x2.*y2-(452+132.*n2+n4).*y4).*ss ...
    ).*fjp - r2.*( (12.*(2+3.*n2).*x6-6.*(68+92.*n2+n4).*x4.*y2+(477+613.*n2+8.*n4).*x2.*y4-(36+44.*n2+n4).*y6).*cs ...
    - 4.*n.*x.*y.*(3.*(34+3.*n2).*x4-(356+37.*n2).*x2.*y2+14.*(8+n2).*y4).*ss ).*fjpp ...
  + r2.*rclose.*(((4.*(3+2.*n2).*x6-2.*(93+50.*n2).*x4.*y2+3.*(69+32.*n2).*x2.*y4-3.*(5+2.*n2).*y6).*cs ...
       - 2.*n.*x.*y.*(2.*(30+n2).*x4-2.*(95+3.*n2).*x2.*y2+(55+2.*n2).*y4).*ss).*fjp3 ...
    - rclose.*(((4+n2).*x6-2.*(25+4.*n2).*x4.*y2+6.*(8+n2).*x2.*y4-3.*y6).*cs - 2.*n.*x.*y.*(11.*x4-28.*x2.*y2+6.*y4).*ss).*fjp4 ...
    + x2.*r2.*(((x4-8.*x2.*y2+6.*y4).*cs - 2.*n.*x.*y.*(x2-2.*y2).*ss).*fjp5 + x2.*y2.*rclose.*cs.*fjp6) ))./r12;
        tmp_nj = reshape(tmp_nj,[],N+1);

        tmp_ni = (-n.*fi.*(n.*(120.*x6-120.*(16+n2).*x4.*y2+(2184+200.*n2+n4).*x2.*y4-20.*(8+n2).*y6).*cs ...
  - 10.*x.*y.*(24.*(3+n2).*x4-2.*(120+56.*n2+n4).*x2.*y2+(72+44.*n2+n4).*y4).*ss) ...
  + rclose.*((24.*(1+4.*n2).*x6-24.*(17+66.*n2+2.*n4).*x4.*y2+(477+1829.*n2+79.*n4).*x2.*y4 ...
       - 4.*(9+34.*n2+2.*n4).*y6).*cs+ 2.*n.*x.*y.*(-24.*(17+3.*n2).*x4+2.*(716+168.*n2+n4).*x2.*y2-(452+132.*n2+n4).*y4).*ss ...
    ).*fip - r2.*( (12.*(2+3.*n2).*x6-6.*(68+92.*n2+n4).*x4.*y2+(477+613.*n2+8.*n4).*x2.*y4-(36+44.*n2+n4).*y6).*cs ...
    - 4.*n.*x.*y.*(3.*(34+3.*n2).*x4-(356+37.*n2).*x2.*y2+14.*(8+n2).*y4).*ss ).*fipp ...
  + r2.*rclose.*(((4.*(3+2.*n2).*x6-2.*(93+50.*n2).*x4.*y2+3.*(69+32.*n2).*x2.*y4-3.*(5+2.*n2).*y6).*cs ...
       - 2.*n.*x.*y.*(2.*(30+n2).*x4-2.*(95+3.*n2).*x2.*y2+(55+2.*n2).*y4).*ss).*fip3 ...
    - rclose.*(((4+n2).*x6-2.*(25+4.*n2).*x4.*y2+6.*(8+n2).*x2.*y4-3.*y6).*cs - 2.*n.*x.*y.*(11.*x4-28.*x2.*y2+6.*y4).*ss).*fip4 ...
    + x2.*r2.*(((x4-8.*x2.*y2+6.*y4).*cs - 2.*n.*x.*y.*(x2-2.*y2).*ss).*fip5 + x2.*y2.*rclose.*cs.*fip6) ))./r12;
        tmp_ni = reshape(tmp_ni,[],N+1);

        sixth_far_xxxxyy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';  

        tmp_nj = (n.*fj.*(n.*x.*y.*(-60.*(14+n2).*x4+(2704+220.*n2+n4).*x2.*y2-60.*(14+n2).*y4).*cs ...
            - 15.*(x2-y2).*(4.*(2+n2).*x4-(112+52.*n2+n4).*x2.*y2+4.*(2+n2).*y4).*ss) ...
  + 3.*rclose.*(x.*y.*(4.*(15+58.*n2+2.*n4).*x4-(195+751.*n2+29.*n4).*x2.*y2+4.*(15+58.*n2+2.*n4).*y4).*cs ...
    + n.*(x2-y2).*(12.*(4+n2).*x4-(668+156.*n2+n4).*x2.*y2+12.*(4+n2).*y4).*ss).*fjp ...
  - 3.*r2.*( x.*y.*((60+80.*n2+n4).*x4-3.*(65+85.*n2+n4).*x2.*y2+(60+80.*n2+n4).*y4).*cs ...
    + n.*(x2-y2).*(3.*(8+n2).*x4-2.*(166+17.*n2).*x2.*y2+3.*(8+n2).*y4).*ss ).*fjpp ...
  + r2.*rclose.*( (3.*x.*y.*((27+14.*n2).*x4-2.*(43+21.*n2).*x2.*y2+(27+14.*n2).*y4).*cs ...
       + n.*(x2-y2).*((20+n2).*x4-(265+8.*n2).*x2.*y2+(20+n2).*y4).*ss).*fjp3 ...
    - 3.*rclose.*((7+n2).*x.*y.*(x4-3.*x2.*y2+y4).*cs+ n.*(x6-14.*x4.*y2+14.*x2.*y4-y6).*ss).*fjp4 ...
    + x.*y.*r2.*( 3.*((x4-3.*x2.*y2+y4).*cs + n.*x.*y.*(-x2+y2).*ss).*fjp5 ...
      + x2.*y2.*rclose.*cs.*fjp6 ) ) )./r12 ;
        tmp_nj = reshape(tmp_nj,[],N+1);

        tmp_ni = (n.*fi.*(n.*x.*y.*(-60.*(14+n2).*x4+(2704+220.*n2+n4).*x2.*y2-60.*(14+n2).*y4).*cs ...
            - 15.*(x2-y2).*(4.*(2+n2).*x4-(112+52.*n2+n4).*x2.*y2+4.*(2+n2).*y4).*ss) ...
  + 3.*rclose.*(x.*y.*(4.*(15+58.*n2+2.*n4).*x4-(195+751.*n2+29.*n4).*x2.*y2+4.*(15+58.*n2+2.*n4).*y4).*cs ...
    + n.*(x2-y2).*(12.*(4+n2).*x4-(668+156.*n2+n4).*x2.*y2+12.*(4+n2).*y4).*ss).*fip ...
  - 3.*r2.*( x.*y.*((60+80.*n2+n4).*x4-3.*(65+85.*n2+n4).*x2.*y2+(60+80.*n2+n4).*y4).*cs ...
    + n.*(x2-y2).*(3.*(8+n2).*x4-2.*(166+17.*n2).*x2.*y2+3.*(8+n2).*y4).*ss ).*fipp ...
  + r2.*rclose.*( (3.*x.*y.*((27+14.*n2).*x4-2.*(43+21.*n2).*x2.*y2+(27+14.*n2).*y4).*cs ...
       + n.*(x2-y2).*((20+n2).*x4-(265+8.*n2).*x2.*y2+(20+n2).*y4).*ss).*fip3 ...
    - 3.*rclose.*((7+n2).*x.*y.*(x4-3.*x2.*y2+y4).*cs+ n.*(x6-14.*x4.*y2+14.*x2.*y4-y6).*ss).*fip4 ...
    + x.*y.*r2.*( 3.*((x4-3.*x2.*y2+y4).*cs + n.*x.*y.*(-x2+y2).*ss).*fip5 ...
      + x2.*y2.*rclose.*cs.*fip6 ) ) )./r12 ;
        tmp_ni = reshape(tmp_ni,[],N+1);

        sixth_far_xxxyyy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';  

        tmp_nj = (-n.*fj.*( n.*(-20.*(8+n2).*x6+(2184+200.*n2+n4).*x4.*y2-120.*(16+n2).*x2.*y4+120.*y6).*cs ...
  + 10.*x.*y.*((72+44.*n2+n4).*x4-2.*(120+56.*n2+n4).*x2.*y2+24.*(3+n2).*y4).*ss ) + rclose.*( ...
      (-4.*(9+34.*n2+2.*n4).*x6+(477+1829.*n2+79.*n4).*x4.*y2-24.*(17+66.*n2+2.*n4).*x2.*y4 ...
       + 24.*(1+4.*n2).*y6).*cs + 2.*n.*x.*y.*((452+132.*n2+n4).*x4-2.*(716+168.*n2+n4).*x2.*y2+24.*(17+3.*n2).*y4).*ss ).*fjp ...
  + r2.*(  ((36+44.*n2+n4).*x6-(477+613.*n2+8.*n4).*x4.*y2+6.*(68+92.*n2+n4).*x2.*y4-12.*(2+3.*n2).*y6).*cs ...
    - 4.*n.*x.*y.*(14.*(8+n2).*x4-(356+37.*n2).*x2.*y2+3.*(34+3.*n2).*y4).*ss ).*fjpp ...
  + r2.*rclose.*(  ((-3.*(5+2.*n2).*x6+3.*(69+32.*n2).*x4.*y2-2.*(93+50.*n2).*x2.*y4+4.*(3+2.*n2).*y6).*cs ...
       + 2.*n.*x.*y.*((55+2.*n2).*x4-2.*(95+3.*n2).*x2.*y2+2.*(30+n2).*y4).*ss).*fjp3 ...
    + rclose.*((3.*x6-6.*(8+n2).*x4.*y2+2.*(25+4.*n2).*x2.*y4-(4+n2).*y6).*cs - 2.*n.*x.*y.*(6.*x4-28.*x2.*y2+11.*y4).*ss).*fjp4 ...
    + y2.*r2.*(((6.*x4-8.*x2.*y2+y4).*cs + 2.*n.*x.*y.*(-2.*x2+y2).*ss).*fjp5 + x2.*y2.*rclose.*cs.*fjp6) ) )./r12;
        tmp_nj = reshape(tmp_nj,[],N+1);

        tmp_ni = (-n.*fi.*( n.*(-20.*(8+n2).*x6+(2184+200.*n2+n4).*x4.*y2-120.*(16+n2).*x2.*y4+120.*y6).*cs ...
  + 10.*x.*y.*((72+44.*n2+n4).*x4-2.*(120+56.*n2+n4).*x2.*y2+24.*(3+n2).*y4).*ss ) + rclose.*( ...
      (-4.*(9+34.*n2+2.*n4).*x6+(477+1829.*n2+79.*n4).*x4.*y2-24.*(17+66.*n2+2.*n4).*x2.*y4 ...
       + 24.*(1+4.*n2).*y6).*cs + 2.*n.*x.*y.*((452+132.*n2+n4).*x4-2.*(716+168.*n2+n4).*x2.*y2+24.*(17+3.*n2).*y4).*ss ).*fip ...
  + r2.*(  ((36+44.*n2+n4).*x6-(477+613.*n2+8.*n4).*x4.*y2+6.*(68+92.*n2+n4).*x2.*y4-12.*(2+3.*n2).*y6).*cs ...
    - 4.*n.*x.*y.*(14.*(8+n2).*x4-(356+37.*n2).*x2.*y2+3.*(34+3.*n2).*y4).*ss ).*fipp ...
  + r2.*rclose.*(  ((-3.*(5+2.*n2).*x6+3.*(69+32.*n2).*x4.*y2-2.*(93+50.*n2).*x2.*y4+4.*(3+2.*n2).*y6).*cs ...
       + 2.*n.*x.*y.*((55+2.*n2).*x4-2.*(95+3.*n2).*x2.*y2+2.*(30+n2).*y4).*ss).*fip3 ...
    + rclose.*((3.*x6-6.*(8+n2).*x4.*y2+2.*(25+4.*n2).*x2.*y4-(4+n2).*y6).*cs - 2.*n.*x.*y.*(6.*x4-28.*x2.*y2+11.*y4).*ss).*fip4 ...
    + y2.*r2.*(((6.*x4-8.*x2.*y2+y4).*cs + 2.*n.*x.*y.*(-2.*x2+y2).*ss).*fip5 + x2.*y2.*rclose.*cs.*fip6) ) )./r12;
        tmp_ni = reshape(tmp_ni,[],N+1);        

        sixth_far_xxyyyy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';  
        
        tmp_nj = (n.*fj.*( n.*x.*y.*((984+140.*n2+n4).*x4-200.*(14+n2).*x2.*y2+600.*y4).*cs ...
            + 5.*((24+20.*n2+n4).*x6-5.*(72+44.*n2+n4).*x4.*y2+120.*(3+n2).*x2.*y4-24.*y6).*ss) ...
            - rclose.*( (5.*x.*y.*((45+169.*n2+11.*n4).*x4-8.*(15+58.*n2+2.*n4).*x2.*y2+24.*(1+4.*n2).*y4).*cs ...
            + n.*((164+60.*n2+n4).*x6-5.*(452+132.*n2+n4).*x4.*y2+120.*(17+3.*n2).*x2.*y4-120.*y6).*ss).*fjp ...
            - 5.*rclose.*(x.*y.*((45+53.*n2+n4).*x4-2.*(60+80.*n2+n4).*x2.*y2+12.*(2+3.*n2).*y4).*cs ...
            + 2.*n.*((8+n2).*x6-14.*(8+n2).*x4.*y2+3.*(34+3.*n2).*x2.*y4-6.*y6).*ss).*fjpp ...
            + r2.*(5.*(2.*x.*y.*(3.*(3+n2).*x4-(27+14.*n2).*x2.*y2+2.*(3+2.*n2).*y4).*cs ...
            + n.*(3.*x6-(55+2.*n2).*x4.*y2+2.*(30+n2).*x2.*y4-4.*y6).*ss).*fjp3 ...
            - y.*( 5.*rclose.*(x.*(3.*x4-2.*(7+n2).*x2.*y2+(4+n2).*y4).*cs - n.*y.*(6.*x4-11.*x2.*y2+y4).*ss).*fjp4 ...
            + y2.*r2.*((5.*x.*(2.*x2-y2).*cs + n.*y.*(-5.*x2+y2).*ss).*fjp5 + x.*y2.*rclose.*cs.*fjp6)) ) ) )./r12;
        tmp_nj = reshape(tmp_nj,[],N+1);        

        tmp_ni = (n.*fi.*( n.*x.*y.*((984+140.*n2+n4).*x4-200.*(14+n2).*x2.*y2+600.*y4).*cs ...
            + 5.*((24+20.*n2+n4).*x6-5.*(72+44.*n2+n4).*x4.*y2+120.*(3+n2).*x2.*y4-24.*y6).*ss) ...
            - rclose.*( (5.*x.*y.*((45+169.*n2+11.*n4).*x4-8.*(15+58.*n2+2.*n4).*x2.*y2+24.*(1+4.*n2).*y4).*cs ...
            + n.*((164+60.*n2+n4).*x6-5.*(452+132.*n2+n4).*x4.*y2+120.*(17+3.*n2).*x2.*y4-120.*y6).*ss).*fip ...
            - 5.*rclose.*(x.*y.*((45+53.*n2+n4).*x4-2.*(60+80.*n2+n4).*x2.*y2+12.*(2+3.*n2).*y4).*cs ...
            + 2.*n.*((8+n2).*x6-14.*(8+n2).*x4.*y2+3.*(34+3.*n2).*x2.*y4-6.*y6).*ss).*fipp ...
            + r2.*(5.*(2.*x.*y.*(3.*(3+n2).*x4-(27+14.*n2).*x2.*y2+2.*(3+2.*n2).*y4).*cs ...
            + n.*(3.*x6-(55+2.*n2).*x4.*y2+2.*(30+n2).*x2.*y4-4.*y6).*ss).*fip3 ...
            - y.*( 5.*rclose.*(x.*(3.*x4-2.*(7+n2).*x2.*y2+(4+n2).*y4).*cs - n.*y.*(6.*x4-11.*x2.*y2+y4).*ss).*fip4 ...
            + y2.*r2.*((5.*x.*(2.*x2-y2).*cs + n.*y.*(-5.*x2+y2).*ss).*fip5 + x.*y2.*rclose.*cs.*fip6)) ) ) )./r12;
        tmp_ni = reshape(tmp_ni,[],N+1);                

        sixth_far_xyyyyy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';  

        tmp_nj = (-n.*x.*fj.*( n.*x.*((184+40.*n2+n4).*x4-300.*(8+n2).*x2.*y2+1800.*y4).*cs ...
             - 30.*y.*((24+20.*n2+n4).*x4-40.*(2+n2).*x2.*y2+24.*y4).*ss ) ...
             + 3.*x.*rclose.*( 5.*x.*((3+11.*n2+n4).*x4-4.*(9+34.*n2+2.*n4).*x2.*y2+24.*(1+4.*n2).*y4).*cs ...
             - 2.*n.*y.*((164+60.*n2+n4).*x4-120.*(4+n2).*x2.*y2+120.*y4).*ss ).*fjp ...
             - 15.*x.*r2.*( x.*(3.*(1+n2).*x4-(36+44.*n2+n4).*x2.*y2+12.*(2+3.*n2).*y4).*cs ...
             - 4.*n.*y.*((8+n2).*x4-3.*(8+n2).*x2.*y2+6.*y4).*ss).*fjpp ...
             + r2.*rclose.*(5.*x.*( 3.*x.*(x4-3.*(5+2.*n2).*x2.*y2+4.*(3+2.*n2).*y4).*cs ...
             + 2.*n.*y.*(-9.*x4+2.*(20+n2).*x2.*y2-12.*y4).*ss ).*fjp3 ...
             + 15.*x.*y2.*rclose.*(x.*(3.*x2-(4+n2).*y2).*cs + 2.*n.*y.*(-2.*x2+y2).*ss ).*fjp4 ...
             + y4.*r2.*(3.*x.*(5.*x.*cs - 2.*n.*y.*ss).*fjp5 + y2.*rclose.*cs.*fjp6) ) )./r12;
        tmp_nj = reshape(tmp_nj,[],N+1);  

        tmp_ni = (-n.*x.*fi.*( n.*x.*((184+40.*n2+n4).*x4-300.*(8+n2).*x2.*y2+1800.*y4).*cs ...
             - 30.*y.*((24+20.*n2+n4).*x4-40.*(2+n2).*x2.*y2+24.*y4).*ss ) ...
             + 3.*x.*rclose.*( 5.*x.*((3+11.*n2+n4).*x4-4.*(9+34.*n2+2.*n4).*x2.*y2+24.*(1+4.*n2).*y4).*cs ...
             - 2.*n.*y.*((164+60.*n2+n4).*x4-120.*(4+n2).*x2.*y2+120.*y4).*ss ).*fip ...
             - 15.*x.*r2.*( x.*(3.*(1+n2).*x4-(36+44.*n2+n4).*x2.*y2+12.*(2+3.*n2).*y4).*cs ...
             - 4.*n.*y.*((8+n2).*x4-3.*(8+n2).*x2.*y2+6.*y4).*ss).*fipp ...
             + r2.*rclose.*(5.*x.*( 3.*x.*(x4-3.*(5+2.*n2).*x2.*y2+4.*(3+2.*n2).*y4).*cs ...
             + 2.*n.*y.*(-9.*x4+2.*(20+n2).*x2.*y2-12.*y4).*ss ).*fip3 ...
             + 15.*x.*y2.*rclose.*(x.*(3.*x2-(4+n2).*y2).*cs + 2.*n.*y.*(-2.*x2+y2).*ss ).*fip4 ...
             + y4.*r2.*(3.*x.*(5.*x.*cs - 2.*n.*y.*ss).*fip5 + y2.*rclose.*cs.*fip6) ) )./r12;
        tmp_ni = reshape(tmp_ni,[],N+1);          

        sixth_far_yyyyyy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';  
  

        sixth_far = cat(3,sixth_far_xxxxxx,sixth_far_xxxxxy,sixth_far_xxxxyy,sixth_far_xxxyyy,sixth_far_xxyyyy,sixth_far_xyyyyy,sixth_far_yyyyyy); 

        sixth(:,iclose,:) = sixth_near + sixth_far;
    end
end

quasi_phase = exp(1i*kappa(:)*nx(:).'*d);


if nargout == 1
    val = quasi_phase.*val;

    if ising == 0 
        for ii = -nsub:nsub
        isub = (abs(nx(:)-ii) > max(ls)) | ifar;

        if any(isub)
        vali = chnk.flex2d.hkdiffgreen(zk,[0;0],[rx(isub).'+ (nx(isub).' - ii)*d;ry(isub).']);
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
        [vali, gradi] = chnk.flex2d.hkdiffgreen(zk,[0;0],[rx(isub).' + (nx(isub).'-ii)*d;ry(isub).']);
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
        [vali, gradi, hessi] = chnk.flex2d.hkdiffgreen(zk,[0;0],[rx(isub).' + (nx(isub).'-ii)*d;ry(isub).']);
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
elseif nargout == 4
    val = quasi_phase.*val;
    grad = quasi_phase.*grad;
    hess = quasi_phase.*hess;
    third = quasi_phase.*third;
    
    if ising == 0
        for ii = -nsub:nsub
        isub = (abs(nx(:)-ii) > max(ls)) | ifar;

        if any(isub)
        [vali, gradi, hessi,thirdi] = chnk.flex2d.hkdiffgreen(zk,[0;0],[rx(isub).' + (nx(isub).'-ii)*d;ry(isub).']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        hessi = reshape(hessi,1,[],3);
        thirdi = reshape(thirdi,1,[],4);

        val(:,isub,:) = val(:,isub,:) - vali.*alpha.^(ii);
        grad(:,isub,:) = grad(:,isub,:) - gradi.*alpha.^(ii);
        hess(:,isub,:) = hess(:,isub,:) - hessi.*alpha.^(ii);
        third(:,isub,:) = third(:,isub,:) - thirdi.*alpha.^(ii);
        end
        end
    end

    val = reshape(val,nkappa*ntarg,nsrc);
    grad = reshape(grad,nkappa*ntarg,nsrc,2);
    hess = reshape(hess,nkappa*ntarg,nsrc,3);
    third = reshape(third,nkappa*ntarg,nsrc,4);
elseif nargout == 5
    val = quasi_phase.*val;
    grad = quasi_phase.*grad;
    hess = quasi_phase.*hess;
    third = quasi_phase.*third;
    fourth = quasi_phase.*fourth;
    
    if ising == 0
        for ii = -nsub:nsub
        isub = (abs(nx(:)-ii) > max(ls)) | ifar;

        if any(isub)
        [vali, gradi, hessi,thirdi,fourthi] = chnk.flex2d.hkdiffgreen(zk,[0;0],[rx(isub).' + (nx(isub).'-ii)*d;ry(isub).']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        hessi = reshape(hessi,1,[],3);
        thirdi = reshape(thirdi,1,[],4);
        fourthi = reshape(fourthi,1,[],5);

        val(:,isub,:) = val(:,isub,:) - vali.*alpha.^(ii);
        grad(:,isub,:) = grad(:,isub,:) - gradi.*alpha.^(ii);
        hess(:,isub,:) = hess(:,isub,:) - hessi.*alpha.^(ii);
        third(:,isub,:) = third(:,isub,:) - thirdi.*alpha.^(ii);
        fourth(:,isub,:) = fourth(:,isub,:) - fourthi.*alpha.^(ii);
        end
        end
    end

    val = reshape(val,nkappa*ntarg,nsrc);
    grad = reshape(grad,nkappa*ntarg,nsrc,2);
    hess = reshape(hess,nkappa*ntarg,nsrc,3);
    third = reshape(third,nkappa*ntarg,nsrc,4);
    fourth = reshape(fourth,nkappa*ntarg,nsrc,5);
elseif nargout == 6
    val = quasi_phase.*val;
    grad = quasi_phase.*grad;
    hess = quasi_phase.*hess;
    third = quasi_phase.*third;
    fourth = quasi_phase.*fourth;
    fifth = quasi_phase.*fifth;
    
    if ising == 0
        for ii = -nsub:nsub
        isub = (abs(nx(:)-ii) > max(ls)) | ifar;

        if any(isub)
        [vali, gradi, hessi,thirdi,fourthi,fifthi] = chnk.flex2d.hkdiffgreen(zk,[0;0],[rx(isub).' + (nx(isub).'-ii)*d;ry(isub).']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        hessi = reshape(hessi,1,[],3);
        thirdi = reshape(thirdi,1,[],4);
        fourthi = reshape(fourthi,1,[],5);
        fifthi = reshape(fifthi,1,[],6);

        val(:,isub,:) = val(:,isub,:) - vali.*alpha.^(ii);
        grad(:,isub,:) = grad(:,isub,:) - gradi.*alpha.^(ii);
        hess(:,isub,:) = hess(:,isub,:) - hessi.*alpha.^(ii);
        third(:,isub,:) = third(:,isub,:) - thirdi.*alpha.^(ii);
        fourth(:,isub,:) = fourth(:,isub,:) - fourthi.*alpha.^(ii);
        fifth(:,isub,:) = fifth(:,isub,:) - fifthi.*alpha.^(ii);
        end
        end
    end

    val = reshape(val,nkappa*ntarg,nsrc);
    grad = reshape(grad,nkappa*ntarg,nsrc,2);
    hess = reshape(hess,nkappa*ntarg,nsrc,3);
    third = reshape(third,nkappa*ntarg,nsrc,4);
    fourth = reshape(fourth,nkappa*ntarg,nsrc,5);
    fifth = reshape(fifth,nkappa*ntarg,nsrc,6);
elseif nargout == 7
    val = quasi_phase.*val;
    grad = quasi_phase.*grad;
    hess = quasi_phase.*hess;
    third = quasi_phase.*third;
    fourth = quasi_phase.*fourth;
    fifth = quasi_phase.*fifth;
    sixth = quasi_phase.*sixth;
    
    if ising == 0
        for ii = -nsub:nsub
        isub = (abs(nx(:)-ii) > max(ls)) | ifar;

        if any(isub)
        [vali, gradi, hessi,thirdi,fourthi,fifthi,sixthi] = chnk.flex2dquas.hkdiffgreen1(zk,[0;0],[rx(isub).' + (nx(isub).'-ii)*d;ry(isub).']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        hessi = reshape(hessi,1,[],3);
        thirdi = reshape(thirdi,1,[],4);
        fourthi = reshape(fourthi,1,[],5);
        fifthi = reshape(fifthi,1,[],6);
        sixthi = reshape(sixthi,1,[],7);

        val(:,isub,:) = val(:,isub,:) - vali.*alpha.^(ii);
        grad(:,isub,:) = grad(:,isub,:) - gradi.*alpha.^(ii);
        hess(:,isub,:) = hess(:,isub,:) - hessi.*alpha.^(ii);
        third(:,isub,:) = third(:,isub,:) - thirdi.*alpha.^(ii);
        fourth(:,isub,:) = fourth(:,isub,:) - fourthi.*alpha.^(ii);
        fifth(:,isub,:) = fifth(:,isub,:) - fifthi.*alpha.^(ii);
        sixth(:,isub,:) = sixth(:,isub,:) - sixthi.*alpha.^(ii);
        end
        end
    end

    val = reshape(val,nkappa*ntarg,nsrc);
    grad = reshape(grad,nkappa*ntarg,nsrc,2);
    hess = reshape(hess,nkappa*ntarg,nsrc,3);
    third = reshape(third,nkappa*ntarg,nsrc,4);
    fourth = reshape(fourth,nkappa*ntarg,nsrc,5);
    fifth = reshape(fifth,nkappa*ntarg,nsrc,6);
    sixth = reshape(sixth,nkappa*ntarg,nsrc,7);
end
end