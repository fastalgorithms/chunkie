function [ds] = helm_axi(rs,drs,dzs,ifun,htables)

    r0s   = sqrt(rs.^2+(rs+drs).^2+dzs.^2);
    alphs = (drs.^2+dzs.^2)./r0s.^2;

    int   = zeros(size(alphs));
    kerns    = int;
    kernsdk  = int;
    kernsda  = int;
    kernsdkk = int;
    kernsdak = int;
    kernsdaa = int;

    iclose = find(alphs<0.2);
    ifar   = 1:numel(alphs);
    ifar(iclose) = [];
    
    alph_far = alphs(ifar);
    r0s_far  = r0s(ifar);
    [i_mid] = find((1-alph_far).*r0s_far < 150);
    imid = ifar(i_mid);
    ifar(i_mid) = [];
    
    alph_mid = alphs(imid);
    r0s_mid  = r0s(imid);
    [i_midnear] = find((1-alph_mid).*r0s_mid < 40);
    imidnear = imid(i_midnear);
    imid(i_midnear) = [];
    
   	alph_midnear = alphs(imidnear);
    r0s_midnear  = r0s(imidnear);
    [i_midnearnear] = find((1-alph_midnear).*r0s_midnear < 40);
    imidnearnear = imidnear(i_midnearnear);
    imidnear(i_midnearnear) = [];
    
    
    [kc,kcdk,kcda,kcdkk,kcdak,kcdaa] ...
        = chnk.axissymhelm2d.helm_axi_close_table(r0s(iclose),alphs(iclose),ifun,htables);
    kerns(iclose) = kc;
    kernsda(iclose) = kcda;
    kernsdk(iclose) = kcdk;
    kernsdkk(iclose) = kcdkk;
    kernsdak(iclose) = kcdak;
    kernsdaa(iclose) = kcdaa;

 	[kf,kfda,kfdk,kfdaa,kfdak,kfdkk] ...
        = chnk.axissymhelm2d.helm_axi_smooth(r0s(imidnearnear),alphs(imidnearnear),...
        ifun,htables.xlege_midnear,htables.wlege_midnear);
    
    kerns(imidnearnear) = kf;
    kernsda(imidnearnear) = kfda;
    kernsdk(imidnearnear) = kfdk;
    kernsdaa(imidnearnear)= kfdaa;
    kernsdak(imidnearnear)= kfdak;
    kernsdkk(imidnearnear)= kfdkk;
    
	[kf,kfda,kfdk,kfdaa,kfdak,kfdkk] ...
        = chnk.axissymhelm2d.helm_axi_smooth(r0s(imidnear),alphs(imidnear),...
        ifun,htables.xlege_midnear,htables.wlege_midnear);
    
    kerns(imidnear) = kf;
    kernsda(imidnear) = kfda;
    kernsdk(imidnear) = kfdk;
    kernsdaa(imidnear)= kfdaa;
    kernsdak(imidnear)= kfdak;
    kernsdkk(imidnear)= kfdkk;
    
    [kf,kfda,kfdk,kfdaa,kfdak,kfdkk] ...
        = chnk.axissymhelm2d.helm_axi_smooth(r0s(imid),alphs(imid),...
        ifun,htables.xlege_mid,htables.wlege_mid);
    
    kerns(imid) = kf;
    kernsda(imid) = kfda;
    kernsdk(imid) = kfdk;
    kernsdaa(imid)= kfdaa;
    kernsdak(imid)= kfdak;
    kernsdkk(imid)= kfdkk;
    
    
    [kf,kfda,kfdk,kfdaa,kfdak,kfdkk] ...
        = chnk.axissymhelm2d.helm_axi_smooth(r0s(ifar),alphs(ifar),...
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
end

