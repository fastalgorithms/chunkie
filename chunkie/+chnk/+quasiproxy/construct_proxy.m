function A = construct_proxy(kh, Proxy,pw,chnkr)

D = kernel('helm', 'd', kh);
S = kernel('helm', 's', kh);

ima=sqrt(-1);

A = (D.eval(Proxy, chnkr)+ima*kh*S.eval(Proxy,chnkr)).*pw.';

return
