function [kerns,kernsda,kernsdk,kernsdaa,kernsdak,kernsdkk, ...
    kerns_i,kernsda_i,kernsdk_i,kernsdaa_i,kernsdak_i,kernsdkk_i] ...
    = helm_axi_smooth_all(r0s,alphs,xlege,wlege)


    [alp,xle] = meshgrid(alphs,xlege);
    xle = cos(xle);
    [r0t,wle] = meshgrid(r0s,wlege);
    rts = sqrt(1-(1-alp).*xle);
    
    efac = exp(1i*r0t.*rts)./rts.*wle;
    r2   = rts.^2;
    k2   = 1i*r0t.*rts;
    
    kerns = sum(chnk.axissymhelm2d.asymint_v(rts,r0t,efac),1);
    kernsdk = sum(chnk.axissymhelm2d.asymintdk_v(rts,r0t,efac),1);
    kernsda = sum(chnk.axissymhelm2d.asymintda_v(xle,rts,r0t,efac),1);
  	kernsdaa = sum(chnk.axissymhelm2d.asymintdaa_v(xle,rts,r0t,efac,k2,r2),1);
 	kernsdak = sum(chnk.axissymhelm2d.asymintdak_v(xle,rts,r0t,efac),1);
  	kernsdkk = sum(chnk.axissymhelm2d.asymintdkk_v(xle,rts,r0t,efac),1);
    
    efac_i = exp(-r0t.*rts)./rts.*wle;
    k2_i = -r0t.*rts;
    r0t_i = r0t*1i;
    
    kerns_i = sum(chnk.axissymhelm2d.asymint_v(rts,r0t_i,efac_i),1);
    kernsdk_i = 1i*sum(chnk.axissymhelm2d.asymintdk_v(rts,r0t_i,efac_i),1);
    kernsda_i = sum(chnk.axissymhelm2d.asymintda_v(xle,rts,r0t_i,efac_i),1);
  	kernsdaa_i = sum(chnk.axissymhelm2d.asymintdaa_v(xle,rts,r0t_i,efac_i,k2_i,r2),1);
 	kernsdak_i = 1i*sum(chnk.axissymhelm2d.asymintdak_v(xle,rts,r0t_i,efac_i),1);
  	kernsdkk_i = -sum(chnk.axissymhelm2d.asymintdkk_v(xle,rts,r0t_i,efac_i),1);


end

