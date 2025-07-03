function uapp = eval_approx(full_sys,chnkr,interface_dens,proxy_dens,numlayer,d,kh,theta,xxtrg)

% will require lots of clean up for multilayer code.


ntot = chnkr.npt;
nproxy = 160;

 sig=interface_dens(1:ntot);
 tau=interface_dens(ntot+1:2*ntot);
 % 
 c1=proxy_dens((numlayer-1)*nproxy+1:(numlayer)*nproxy);


%  create the approximate solution at the target point
D = kernel('helm', 'd', kh);
S = kernel('helm', 's', kh);


% make neighbors
chnkr_l = chnkr;
chnkr_l.r(1,:) = chnkr_l.r(1,:)-d;
chnkr_r = chnkr;
chnkr_r.r(1,:) = chnkr_r.r(1,:)+d;

ww = chnkr.wts;
ww = reshape(ww,1,numel(ww));
wwL = chnkr_l.wts;
wwL = reshape(wwL,1,numel(wwL));
wwR = chnkr_r.wts;
wwR = reshape(wwR,1,numel(wwR));

DL = D.eval(chnkr_l,xxtrg).*wwL;
DD = D.eval(chnkr,xxtrg).*ww;
DR = D.eval(chnkr_r,xxtrg).*wwR;

SL = S.eval(chnkr_l,xxtrg).*wwL;
SS  = S.eval(chnkr,xxtrg).*ww;
SR = S.eval(chnkr_r,xxtrg).*wwR;

% bloch phase
ima = sqrt(-1);
alpha = exp(ima*kh*d*cos(theta));


us = alpha^(-1)*SL*tau+alpha*SR*tau+SS*tau;

ud = alpha^(-1)*DL*sig+alpha*DR*sig+DD*sig;


uproxy = chnk.quasiproxy.construct_proxy(kh,full_sys.Cproxy{numlayer},full_sys.Cproxy{end},xxtrg)*c1;


uapp = us+ud+uproxy;

return


