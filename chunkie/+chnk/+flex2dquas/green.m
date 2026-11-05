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
%          CHNK.FLEX2D.HKDIFFGREEN
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


nx = round(rx/d);
rx = rx - nx*d;


rx2 = rx.*rx;
ry2 = ry.*ry;

r2 = rx2+ry2;

r = sqrt(r2);

npt = size(r,1);

ythresh = d;
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
        [vali,gradi,hessi,thirdi,fourthi,fifthi,sixthi] = chnk.flex2d.hkdiffgreen(zk,[0;0],[rxi.';ryclose.']);
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

    nord = max(nargout,1)-1;
    Dj = expansion_derivs(sn(:,:,1),zk,rxclose,ryclose,rclose,nord);
    Di = expansion_derivs(sn(:,:,2),1i*zk,rxclose,ryclose,rclose,nord);

    val(:,iclose) = val_near + Dj{1} - Di{1};
    if nargout > 1
        grad(:,iclose,:) = grad_near + Dj{2} - Di{2};
    end
    if nargout > 2
        hess(:,iclose,:) = hess_near + Dj{3} - Di{3};
    end
    if nargout > 3
        third(:,iclose,:) = third_near + Dj{4} - Di{4};
    end
    if nargout > 4
        fourth(:,iclose,:) = fourth_near + Dj{5} - Di{5};
    end
    if nargout > 5
        fifth(:,iclose,:) = fifth_near + Dj{6} - Di{6};
    end
    if nargout > 6
        sixth(:,iclose,:) = sixth_near + Dj{7} - Di{7};
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
        [vali, gradi, hessi,thirdi,fourthi,fifthi,sixthi] = chnk.flex2d.hkdiffgreen(zk,[0;0],[rx(isub).' + (nx(isub).'-ii)*d;ry(isub).']);
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

function D = expansion_derivs(c,k,rx,ry,r,nord)
% derivatives of orders 0...nord of the local expansion
%   f = 0.25i c_0 J_0(k r) + 0.5i sum_{n>=1} c_n J_n(k r) cos(n theta)
% D{o+1} is (nkappa, npt, o+1) with components ordered x...x, ..., y...y.
%
% Writes f = sum_j b_j u_j with u_j = J_j(k r) e^{i j theta}, and uses
% (d_x + i d_y) u_j = -k u_{j+1} and (d_x - i d_y) u_j = k u_{j-1}, which
% involve no division by r.

N = size(c,2)-1;
jmax = N + nord;
npt = numel(r);

eip = (rx(:)+1i*ry(:))./r(:);
eip(r(:) < 1e-14) = 1;
Jm = zeros(npt,jmax+1);
for m = 0:jmax
    Jm(:,m+1) = besselj(m,k*r(:));
end
js = -jmax:jmax;
u = Jm(:,abs(js)+1).*((-1).^(js.*(js<0))).*eip.^js;   % (npt, 2*jmax+1)

ns = -N:N;
b = 0.25i*c(:,abs(ns)+1).*((-1).^(ns.*(ns<0)));         % (nkappa, 2N+1)

D = cell(nord+1,1);
for o = 0:nord
    % U(:,:,p+1) = sum_j b_j u_{j+2p-o}
    U = zeros(size(c,1),npt,o+1);
    for p = 0:o
        U(:,:,p+1) = b*u(:,ns+2*p-o+jmax+1).';
    end
    Do = zeros(size(c,1),npt,o+1);
    for t = 0:o
        a = o-t;
        for p = 0:o
            % coefficient of D+^p D-^(o-p) in d_x^a d_y^t
            w = 0;
            for sx = max(0,p-t):min(a,p)
                ty = p-sx;
                w = w + nchoosek(a,sx)*nchoosek(t,ty)*(-1)^(t-ty);
            end
            w = w/(2^a*(2i)^t);
            if w ~= 0
                Do(:,:,t+1) = Do(:,:,t+1) + w*(-1)^p*k^o*U(:,:,p+1);
            end
        end
    end
    D{o+1} = Do;
end
end
