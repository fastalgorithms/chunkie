function [kerns,kernsda,kernsdk,kernsdaa,kernsdak,kernsdkk] ...
    = helm_axi_smooth2(r0s,alphs,ifun,xlege,wlege)

 	int   = zeros(size(alphs));
    kerns = int;
    kernsdk = int;
    kernsda = int;
    kernsdkk= int;
    kernsdak= int;
    kernsdaa= int;
       
    for ii=1:numel(alphs)
        
        alph = alphs(ii);
        r0   = r0s(ii);
        
        if (ifun == 1)
            rt = sqrt(1-(1-alph)*xlege);
            kern = sum(chnk.axissymhelm2d.asymint(xlege,alph,r0).*wlege);
            kerndk = sum(chnk.axissymhelm2d.asymintdk(xlege,alph,r0).*wlege);
            kernda = sum(chnk.axissymhelm2d.asymintda(xlege,alph,r0).*wlege);
            kerndak= sum(chnk.axissymhelm2d.asymintdak(xlege,alph,r0).*wlege);
            kerndaa= sum(chnk.axissymhelm2d.asymintdaa(xlege,alph,r0).*wlege);
            kerndkk= sum(chnk.axissymhelm2d.asymintdkk(xlege,alph,r0).*wlege);
        end
        if (ifun == 2)
            rt = sqrt(1-(1-alph)*xlege);
            kern = sum(chnk.axissymhelm2d.asymint(xlege,alph,1i*r0).*wlege);
            kerndk = 1i*sum(chnk.axissymhelm2d.asymintdk(xlege,alph,1i*r0).*wlege);
            kernda = sum(chnk.axissymhelm2d.asymintda(xlege,alph,1i*r0).*wlege);
            kerndak= 1i*sum(chnk.axissymhelm2d.asymintdak(xlege,alph,1i*r0).*wlege);
            kerndaa= sum(chnk.axissymhelm2d.asymintdaa(xlege,alph,1i*r0).*wlege);
            kerndkk= -sum(chnk.axissymhelm2d.asymintdkk(xlege,alph,1i*r0).*wlege);
        end
        if (ifun == 3)
            rt = sqrt(1-(1-alph)*xlege);
            kern = sum(chnk.axissymhelm2d.asymint(xlege,alph,r0).*wlege);
            kerndk = sum(chnk.axissymhelm2d.asymintdk(xlege,alph,r0).*wlege);
            kernda = sum(chnk.axissymhelm2d.asymintda(xlege,alph,r0).*wlege);
            kerndak= sum(chnk.axissymhelm2d.asymintdak(xlege,alph,r0).*wlege);
            kerndaa= sum(chnk.axissymhelm2d.asymintdaa(xlege,alph,r0).*wlege);
            kerndkk= sum(chnk.axissymhelm2d.asymintdkk(xlege,alph,r0).*wlege);
            
            kern2 = sum(chnk.axissymhelm2d.asymint(xlege,alph,1i*r0).*wlege);
            kern2dk = 1i*sum(chnk.axissymhelm2d.asymintdk(xlege,alph,1i*r0).*wlege);
            kern2da = sum(chnk.axissymhelm2d.asymintda(xlege,alph,1i*r0).*wlege);
            kern2dak= 1i*sum(chnk.axissymhelm2d.asymintdak(xlege,alph,1i*r0).*wlege);
            kern2daa= sum(chnk.axissymhelm2d.asymintdaa(xlege,alph,1i*r0).*wlege);
            kern2dkk= -sum(chnk.axissymhelm2d.asymintdkk(xlege,alph,1i*r0).*wlege);
            
            kern    = kern    - kern2;
            kernda  = kernda  - kern2da;
            kerndk  = kerndk  - kern2dk;
            kerndaa = kerndaa - kern2daa;
            kerndak = kerndak - kern2dak;
            kerndkk = kerndkk - kern2dkk;
        end
        
        kerns(ii) = kern;
        kernsdk(ii) = kerndk;
        kernsda(ii) = kernda;
        kernsdaa(ii)= kerndaa;
        kernsdak(ii)= kerndak;
        kernsdkk(ii)= kerndkk;
        
    end
    
    
    
end

