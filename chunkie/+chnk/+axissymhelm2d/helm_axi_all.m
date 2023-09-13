function [dsk, dsik, dsdiff] = helm_axi_all(rs,drs,dzs,htables)

    r0s   = sqrt(rs.^2+(rs+drs).^2+dzs.^2);
    alphs = (drs.^2+dzs.^2)./r0s.^2;

    int   = zeros(size(alphs));
    kerns_k    = int;
    kernsdk_k  = int;
    kernsda_k  = int;
    kernsdkk_k = int;
    kernsdak_k = int;
    kernsdaa_k = int;

    kerns_ik    = int;
    kernsdk_ik  = int;
    kernsda_ik  = int;
    kernsdkk_ik = int;
    kernsdak_ik = int;
    kernsdaa_ik = int;
    
    
    kerns_diff    = int;
    kernsdk_diff  = int;
    kernsda_diff  = int;
    kernsdkk_diff = int;
    kernsdak_diff = int;
    kernsdaa_diff = int;

    
    
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
        = chnk.axissymhelm2d.helm_axi_close_table(r0s(iclose),alphs(iclose),1,htables);
    kerns_k(iclose) = kc;
    kernsda_k(iclose) = kcda;
    kernsdk_k(iclose) = kcdk;
    kernsdkk_k(iclose) = kcdkk;
    kernsdak_k(iclose) = kcdak;
    kernsdaa_k(iclose) = kcdaa;

        
    [kc,kcdk,kcda,kcdkk,kcdak,kcdaa] ...
        = chnk.axissymhelm2d.helm_axi_close_table(r0s(iclose),alphs(iclose),2,htables);
    kerns_ik(iclose) = kc;
    kernsda_ik(iclose) = kcda;
    kernsdk_ik(iclose) = kcdk;
    kernsdkk_ik(iclose) = kcdkk;
    kernsdak_ik(iclose) = kcdak;
    kernsdaa_ik(iclose) = kcdaa;

            
    [kc,kcdk,kcda,kcdkk,kcdak,kcdaa] ...
        = chnk.axissymhelm2d.helm_axi_close_table(r0s(iclose),alphs(iclose),3,htables);
    kerns_diff(iclose) = kc;
    kernsda_diff(iclose) = kcda;
    kernsdk_diff(iclose) = kcdk;
    kernsdkk_diff(iclose) = kcdkk;
    kernsdak_diff(iclose) = kcdak;
    kernsdaa_diff(iclose) = kcdaa;

    
    
 	[kf_k,kfda_k,kfdk_k,kfdaa_k,kfdak_k,kfdkk_k, ...
     kf_ik,kfda_ik,kfdk_ik,kfdaa_ik,kfdak_ik,kfdkk_ik] ...
        = chnk.axissymhelm2d.helm_axi_smooth_all(r0s(imidnearnear),alphs(imidnearnear),...
        htables.xlege_midnear,htables.wlege_midnear);
    
    kerns_k(imidnearnear) = kf_k;
    kernsda_k(imidnearnear) = kfda_k;
    kernsdk_k(imidnearnear) = kfdk_k;
    kernsdaa_k(imidnearnear)= kfdaa_k;
    kernsdak_k(imidnearnear)= kfdak_k;
    kernsdkk_k(imidnearnear)= kfdkk_k;

        
    kerns_ik(imidnearnear) = kf_ik;
    kernsda_ik(imidnearnear) = kfda_ik;
    kernsdk_ik(imidnearnear) = kfdk_ik;
    kernsdaa_ik(imidnearnear)= kfdaa_ik;
    kernsdak_ik(imidnearnear)= kfdak_ik;
    kernsdkk_ik(imidnearnear)= kfdkk_ik;

            
    kerns_diff(imidnearnear) = kf_k - kf_ik;
    kernsda_diff(imidnearnear) = kfda_k - kfda_ik;
    kernsdk_diff(imidnearnear) = kfdk_k - kfdk_ik;
    kernsdaa_diff(imidnearnear)= kfdaa_k - kfdaa_ik;
    kernsdak_diff(imidnearnear)= kfdak_k - kfdak_ik;
    kernsdkk_diff(imidnearnear)= kfdkk_k - kfdkk_ik;

    
	[kf_k,kfda_k,kfdk_k,kfdaa_k,kfdak_k,kfdkk_k, ...
     kf_ik,kfda_ik,kfdk_ik,kfdaa_ik,kfdak_ik,kfdkk_ik] ...
        = chnk.axissymhelm2d.helm_axi_smooth_all(r0s(imidnear),alphs(imidnear),...
        htables.xlege_midnear,htables.wlege_midnear);
    
    kerns_k(imidnear) = kf_k;
    kernsda_k(imidnear) = kfda_k;
    kernsdk_k(imidnear) = kfdk_k;
    kernsdaa_k(imidnear)= kfdaa_k;
    kernsdak_k(imidnear)= kfdak_k;
    kernsdkk_k(imidnear)= kfdkk_k;
    
            
    kerns_ik(imidnear) = kf_ik;
    kernsda_ik(imidnear) = kfda_ik;
    kernsdk_ik(imidnear) = kfdk_ik;
    kernsdaa_ik(imidnear)= kfdaa_ik;
    kernsdak_ik(imidnear)= kfdak_ik;
    kernsdkk_ik(imidnear)= kfdkk_ik;

            
    kerns_diff(imidnear) = kf_k - kf_ik;
    kernsda_diff(imidnear) = kfda_k - kfda_ik;
    kernsdk_diff(imidnear) = kfdk_k - kfdk_ik;
    kernsdaa_diff(imidnear)= kfdaa_k - kfdaa_ik;
    kernsdak_diff(imidnear)= kfdak_k - kfdak_ik;
    kernsdkk_diff(imidnear)= kfdkk_k - kfdkk_ik;

    
    [kf_k,kfda_k,kfdk_k,kfdaa_k,kfdak_k,kfdkk_k, ...
     kf_ik,kfda_ik,kfdk_ik,kfdaa_ik,kfdak_ik,kfdkk_ik] ...
        = chnk.axissymhelm2d.helm_axi_smooth_all(r0s(imid),alphs(imid),...
        htables.xlege_mid,htables.wlege_mid);
    
    kerns_k(imid) = kf_k;
    kernsda_k(imid) = kfda_k;
    kernsdk_k(imid) = kfdk_k;
    kernsdaa_k(imid)= kfdaa_k;
    kernsdak_k(imid)= kfdak_k;
    kernsdkk_k(imid)= kfdkk_k;
              
    kerns_ik(imid) = kf_ik;
    kernsda_ik(imid) = kfda_ik;
    kernsdk_ik(imid) = kfdk_ik;
    kernsdaa_ik(imid)= kfdaa_ik;
    kernsdak_ik(imid)= kfdak_ik;
    kernsdkk_ik(imid)= kfdkk_ik;

            
    kerns_diff(imid) = kf_k - kf_ik;
    kernsda_diff(imid) = kfda_k - kfda_ik;
    kernsdk_diff(imid) = kfdk_k - kfdk_ik;
    kernsdaa_diff(imid)= kfdaa_k - kfdaa_ik;
    kernsdak_diff(imid)= kfdak_k - kfdak_ik;
    kernsdkk_diff(imid)= kfdkk_k - kfdkk_ik;

    
   [kf_k,kfda_k,kfdk_k,kfdaa_k,kfdak_k,kfdkk_k, ...
     kf_ik,kfda_ik,kfdk_ik,kfdaa_ik,kfdak_ik,kfdkk_ik] ...
        = chnk.axissymhelm2d.helm_axi_smooth_all(r0s(ifar),alphs(ifar),...
        htables.xlege,htables.wlege);
    
    kerns_k(ifar) = kf_k;
    kernsda_k(ifar) = kfda_k;
    kernsdk_k(ifar) = kfdk_k;
    kernsdaa_k(ifar)= kfdaa_k;
    kernsdak_k(ifar)= kfdak_k;
    kernsdkk_k(ifar)= kfdkk_k;
    
                  
    kerns_ik(ifar) = kf_ik;
    kernsda_ik(ifar) = kfda_ik;
    kernsdk_ik(ifar) = kfdk_ik;
    kernsdaa_ik(ifar)= kfdaa_ik;
    kernsdak_ik(ifar)= kfdak_ik;
    kernsdkk_ik(ifar)= kfdkk_ik;

            
    kerns_diff(ifar) = kf_k - kf_ik;
    kernsda_diff(ifar) = kfda_k - kfda_ik;
    kernsdk_diff(ifar) = kfdk_k - kfdk_ik;
    kernsdaa_diff(ifar)= kfdaa_k - kfdaa_ik;
    kernsdak_diff(ifar)= kfdak_k - kfdak_ik;
    kernsdkk_diff(ifar)= kfdkk_k - kfdkk_ik;


    [doutk] = chnk.axissymhelm2d.der_ak_2_grad(rs,drs,dzs,kerns_k,kernsda_k,kernsdk_k,...
                 kernsdaa_k,kernsdak_k,kernsdkk_k);
    [dsk] = chnk.axissymhelm2d.div_by_kap(rs,drs,dzs,doutk);
    dsk.intdrz = -dsk.intdrz;
     dsk.intdzz = -dsk.intdzz;
     
    
    [doutik] = chnk.axissymhelm2d.der_ak_2_grad(rs,drs,dzs,kerns_ik,kernsda_ik,kernsdk_ik,...
                 kernsdaa_ik,kernsdak_ik,kernsdkk_ik);
    [dsik] = chnk.axissymhelm2d.div_by_kap(rs,drs,dzs,doutik);
    dsik.intdrz = -dsik.intdrz;
     dsik.intdzz = -dsik.intdzz; 
     
     
         
    [doutdiff] = chnk.axissymhelm2d.der_ak_2_grad(rs,drs,dzs,kerns_diff,kernsda_diff,kernsdk_diff,...
                 kernsdaa_diff,kernsdak_diff,kernsdkk_diff);
    [dsdiff] = chnk.axissymhelm2d.div_by_kap(rs,drs,dzs,doutdiff);
    dsdiff.intdrz = -dsdiff.intdrz;
     dsdiff.intdzz = -dsdiff.intdzz; 

end

