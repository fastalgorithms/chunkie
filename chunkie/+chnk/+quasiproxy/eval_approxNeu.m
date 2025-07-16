function uapp = eval_approxNeu(full_sys,chnkr,interface_dens,proxy_dens,bragg_coef,KK,ht,d,theta,xxtrg)

% will require lots of clean up for multilayer code.

% detect layer

numlayer = 1;

kh = full_sys.khs(numlayer);

nproxy = size(full_sys.Cproxy{numlayer}.r,2);

if xxtrg.r(2,1)>ht
    % use bragg expansion
    topcoef = bragg_coef(1:KK);
    kh1 = full_sys.khs(1);


    uapp = 0;
    n = (KK-1)/2;
    ima = sqrt(-1);
    for j= 1:KK
        kappa = kh1*cos(theta)+2*pi*(j-n);
        knu = sqrt(kh1^2-kappa^2);    
        uapp = uapp+topcoef(j)*exp(ima*kappa*xxtrg.r(1))*exp(ima*knu*(xxtrg.r(2)-ht));
    end

    
else    
% use integral operator form of solution    
% extract the different parts of the solution to the linear system.
 c1=proxy_dens(1:nproxy);


 %  create the approximate solution at the target point

uproxy = chnk.quasiproxy.construct_proxy(kh,full_sys.Cproxy{numlayer},full_sys.Cproxy{end},xxtrg)*c1;
uproxy = uproxy(1,:);

% make neighbors
chnkr_l = chnkr;
chnkr_l.r(1,:) = chnkr_l.r(1,:)-d;
chnkr_r = chnkr;
chnkr_r.r(1,:) = chnkr_r.r(1,:)+d;

C = @(s,t) chnk.helm2d.kern(kh,s,t,'s');
% get weights
ww = chnkr.wts; ww = ww(:).'; 
% duplicate weights for second density
ww = repmat(ww,2,1); ww = ww(:).';

% evaluate layer potentials
TrL = C(chnkr_l,xxtrg).*ww;
TrM = C(chnkr  ,xxtrg).*ww;
TrR = C(chnkr_r,xxtrg).*ww;

% bloch phase
ima = sqrt(-1);
alpha = exp(ima*kh*d*cos(theta));
udens = alpha^(-1)*(TrL*interface_dens) + (TrM*interface_dens) ...
    + alpha*(TrR*interface_dens);

uapp = udens+uproxy;
end


return
end






function ireg = find_layer(cgrph_lab,hb,ht,xxtrg)

ireg = chunkgraphinregion(cgrph_lab,xxtrg);

% ireg(ireg==2) = 1;
% ireg(ireg==1) = 1+(pt(2,ireg==1).' < hb); 
% ireg(ireg==3) = 2;

ireg = ireg - 1;
ireg(ireg == 0) = 3 + (xxtrg(2,ireg==0) < hb);



end