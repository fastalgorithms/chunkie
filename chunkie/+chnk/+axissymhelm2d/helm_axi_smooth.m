function [kerns,kernsda,kernsdk,kernsdaa,kernsdak,kernsdkk,...
    varargout] ...
    = helm_axi_smooth(r0s,alphs,ifun,xlege,wlege)


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

        if (nargout == 12)
            varargout{1} = kern2_v;
            varargout{2} = kern2da_v;
            varargout{3} = 1i*kern2dk_v;
            varargout{4} = kern2daa_v;
            varargout{5} = 1i*kern2dak_v;
            varargout{6} = -kern2dkk_v;
        end
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





