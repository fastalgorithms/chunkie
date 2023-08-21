function [kerns,kernsda,kernsdk,kernsdaa,kernsdak,kernsdkk] ...
    = helm_axi_smooth2(r0s,alphs,ifun,xlege,wlege)


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
    
    kerns = sum(chnk.axissymhelm2d.asymint_v(rts,r0t,efac),1);
    kernsdk = sum(chnk.axissymhelm2d.asymintdk_v(rts,r0t,efac),1);
    kernsda = sum(chnk.axissymhelm2d.asymintda_v(xle,rts,r0t,efac),1);
  	kernsdaa = sum(chnk.axissymhelm2d.asymintdaa_v(xle,rts,r0t,efac,k2,r2),1);
 	kernsdak = sum(chnk.axissymhelm2d.asymintdak_v(xle,rts,r0t,efac),1);
  	kernsdkk = sum(chnk.axissymhelm2d.asymintdkk_v(xle,rts,r0t,efac),1);
    
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
        kern2_v = sum(chnk.axissymhelm2d.asymint_v(rts,r0t,efac),1);
        kern2dk_v = sum(chnk.axissymhelm2d.asymintdk_v(rts,r0t,efac),1);
        kern2da_v = sum(chnk.axissymhelm2d.asymintda_v(xle,rts,r0t,efac),1);
        kern2daa_v = sum(chnk.axissymhelm2d.asymintdaa_v(xle,rts,r0t,efac,k2,r2),1);
        kern2dak_v = sum(chnk.axissymhelm2d.asymintdak_v(xle,rts,r0t,efac),1);
        kern2dkk_v = sum(chnk.axissymhelm2d.asymintdkk_v(xle,rts,r0t,efac),1);
        kerns      = kerns       -    kern2_v;
        kernsdk    = kernsdk     - 1i*kern2dk_v;
        kernsda    = kernsda     -    kern2da_v;
        kernsdak   = kernsdak    - 1i*kern2dak_v;
        kernsdaa   = kernsdaa    -    kern2daa_v;
        kernsdkk   = kernsdkk    +    kern2dkk_v;     
    end

end

