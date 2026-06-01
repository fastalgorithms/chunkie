function [val,vald] = psi_eval(x,dsig)
    z = x.*x./(2*dsig.*dsig);
    vexp = chnk.smoother.expeval(z,1);
    val  = vexp/(4*pi);
    vald = -exp(-z)./(2*pi*x);
end
