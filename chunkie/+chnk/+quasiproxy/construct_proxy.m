function A = construct_proxy(kh, Proxy,pw,chnkr)

repcoef = [1,1i*kh];
A = chnk.helm2d.kern(kh,Proxy,chnkr,'c',repcoef).*pw.';

return
