function [sout] ...
    = helm_axi_smooth(r0s,alphs,ifun,xlege,wlege)


%%%%
%
%       On output, a cell array with entries
%       *kerns
%       *kernsda
%       *kernsdk
%       *kernsdkk
%       *kernsdak
%       *kernsdaa
%
%%%%


    [alp,xle] = meshgrid(alphs,xlege);
    xle = cos(xle);
    [r0t,wle] = meshgrid(r0s,wlege);
    rts = sqrt(1-(1-alp).*xle);
    
    if (ifun == 2)
        r0t = r0t*1i;
    end
    
    efac = exp(1i*r0t.*rts)./rts.*wle;
    r2   = rts.^2;
    k2   = 1i*r0t.*rts;
    
    kerns = sum(asymint_v(rts,r0t,efac),1);
    kernsdk = sum(asymintdk_v(rts,r0t,efac),1);
    kernsda = sum(asymintda_v(xle,rts,r0t,efac),1);
  	kernsdaa = sum(asymintdaa_v(xle,rts,r0t,efac,k2,r2),1);
 	kernsdak = sum(asymintdak_v(xle,rts,r0t,efac),1);
  	kernsdkk = sum(asymintdkk_v(xle,rts,r0t,efac),1);
    
    if (ifun == 2)
        kernsdk =   kernsdk*1i;
        kernsdak =  kernsdak*1i;
        kernsdkk = -kernsdkk;
    end
    
    if (ifun == 3)
        r0t = r0t*1i;
       	efac = exp(1i*r0t.*rts)./rts.*wle;
        r2   = rts.^2;
        k2   = 1i*r0t.*rts;
        kern2_v = sum(asymint_v(rts,r0t,efac),1);
        kern2dk_v = sum(asymintdk_v(rts,r0t,efac),1);
        kern2da_v = sum(asymintda_v(xle,rts,r0t,efac),1);
        kern2daa_v = sum(asymintdaa_v(xle,rts,r0t,efac,k2,r2),1);
        kern2dak_v = sum(asymintdak_v(xle,rts,r0t,efac),1);
        kern2dkk_v = sum(asymintdkk_v(xle,rts,r0t,efac),1);
        kerns      = kerns       -    kern2_v;
        kernsdk    = kernsdk     - 1i*kern2dk_v;
        kernsda    = kernsda     -    kern2da_v;
        kernsdak   = kernsdak    - 1i*kern2dak_v;
        kernsdaa   = kernsdaa    -    kern2daa_v;
        kernsdkk   = kernsdkk    +    kern2dkk_v;     
    end

    s ={};
    s{1}    = kerns;
    s{3}  = kernsdk;
    s{2}  = kernsda;
    s{5} = kernsdak;
    s{6} = kernsdaa;
    s{4} = kernsdkk;

    if (ifun == 1 || ifun==4)
        sout.rk = s;
    end
    if (ifun ==2)
        sout.ik = s;
    end
    if (ifun ==3)
        sout.dk = s;
    end
        
    if (ifun == 4)
        r0t = r0t*1i;
        efac = exp(1i*r0t.*rts)./rts.*wle;
        r2   = rts.^2;
        k2   = 1i*r0t.*rts;
        kern2_v = sum(asymint_v(rts,r0t,efac),1);
        kern2dk_v = sum(asymintdk_v(rts,r0t,efac),1);
        kern2da_v = sum(asymintda_v(xle,rts,r0t,efac),1);
        kern2daa_v = sum(asymintdaa_v(xle,rts,r0t,efac,k2,r2),1);
        kern2dak_v = sum(asymintdak_v(xle,rts,r0t,efac),1);
        kern2dkk_v = sum(asymintdkk_v(xle,rts,r0t,efac),1);

        ik = {};
        ik{1} = kern2_v;
        ik{2} = kern2da_v;
        ik{3} = 1i*kern2dk_v;
        ik{6} = kern2daa_v;
        ik{5} = 1i*kern2dak_v;
        ik{4} = -kern2dkk_v;
        sout.ik = ik;

        dk = {};
        for ii=1:6
            dk{ii} = sout.rk{ii}-sout.ik{ii};
        end
        sout.dk = dk;
    end

end

function [val] = asymint_v(r,k,efac)
    val = efac;
end

function [val] = asymintda_v(x,r,k,efac)
    val = x.*efac./(r.^2).*(1i*k.*r-1)/2;
    %val = exp(1i*k*sqrt(1-(1-a)*cos(x)))./sqrt(1-(1-a)*cos(x));
end

function [val] = asymintdaa_v(x,r,k,efac,k2,r2)
    rdiv = r2.^2;
    prefac = 0.25*(x.^2).*efac./rdiv;
    fac = k2.^2-3*k2+3;
    val = fac.*prefac;
    %val = -exp(1i*k*sqrt(1-(1-a)*cos(x))).*sqrt(1-(1-a)*cos(x));
end

function [val] = asymintdak_v(x,r,k,efac)
    val = -(k/2).*x.*efac;
end



function [val] = asymintdk_v(r,k,efac)
    val = 1i*efac.*r;
end


function [val] = asymintdkk_v(x,r,k,efac)
    val = -efac.*(r.^2);
end





