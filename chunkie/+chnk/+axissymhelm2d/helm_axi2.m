function [ds] = helm_axi2(rs,drs,dzs,ifun,htables)

    r0s   = sqrt(rs.^2+(rs+drs).^2+dzs.^2);
    alphs = (drs.^2+dzs.^2)./r0s.^2;

    int   = zeros(size(alphs));
    kerns = int;
    kernsdk = int;
    kernsda = int;
    kernsdkk= int;
    kernsdak= int;
    kernsdaa= int;

    iclose = find(alphs<0.2);
    ifar   = 1:numel(alphs);
    ifar(iclose) = [];
    
    alph_far = alphs(ifar);
    r0s_far  = r0s(ifar);
    [i_mid] = find(alph_far.*r0s_far < 150);
    imid = ifar(i_mid);
    ifar(i_mid) = [];
    
    [kc,kcdk,kcda,kcdkk,kcdak,kcdaa] ...
        = helm_axi_close_table(r0s(iclose),alphs(iclose),ifun,htables);
    kerns(iclose) = kc;
    kernsda(iclose) = kcda;
    kernsdk(iclose) = kcdk;
    kernsdkk(iclose) = kcdkk;
    kernsdak(iclose) = kcdak;
    kernsdaa(iclose) = kcdaa;
 
    [kf,kfda,kfdk,kfdaa,kfdak,kfdkk] ...
        = chnk.axissymhelm2d.helm_axi_smooth2(r0s(imid),alphs(imid),...
        ifun,htables.xlege_mid,htables.wlege_mid);
    
    kerns(imid) = kf;
    kernsda(imid) = kfda;
    kernsdk(imid) = kfdk;
    kernsdaa(imid)= kfdaa;
    kernsdak(imid)= kfdak;
    kernsdkk(imid)= kfdkk;
    
    
    [kf,kfda,kfdk,kfdaa,kfdak,kfdkk] ...
        = chnk.axissymhelm2d.helm_axi_smooth2(r0s(ifar),alphs(ifar),...
        ifun,htables.xlege,htables.wlege);
    
    kerns(ifar) = kf;
    kernsda(ifar) = kfda;
    kernsdk(ifar) = kfdk;
    kernsdaa(ifar)= kfdaa;
    kernsdak(ifar)= kfdak;
    kernsdkk(ifar)= kfdkk;

    [dout] = chnk.axissymhelm2d.der_ak_2_grad(rs,drs,dzs,kerns,kernsda,kernsdk,...
                 kernsdaa,kernsdak,kernsdkk);
    [ds] = chnk.axissymhelm2d.div_by_kap(rs,drs,dzs,dout);
    ds.intdrz = -ds.intdrz;
     ds.intdzz = -ds.intdzz;
     
%     int(ii)   = ds.int;
%     intdq(ii) = ds.intdq;
%     intdr(ii) = ds.intdr;
%     intdz(ii) = ds.intdz;
%     intdzz(ii)= ds.intdzz;
%     intdrq(ii)= ds.intdrq;
%     intdrz(ii)= ds.intdrz;
%     intdqz(ii)= ds.intdqz;
%     
%     ds = [];
%     ds.int = int;
%     ds.intdq = intdq;
%     ds.intdz = intdz;
%     ds.intdr = intdr;
%     ds.intdzz = intdzz;
%     ds.intdrq = intdrq;
%     ds.intdrz = intdrz;
%     ds.intdqz = intdqz;
end

