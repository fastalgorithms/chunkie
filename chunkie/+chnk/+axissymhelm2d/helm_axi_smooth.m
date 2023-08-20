function [kern,kernda,kerndk,kerndaa,kerndak,kerndkk] ...
    = helm_axi_smooth(r0,alph,ifun,xlege,wlege)

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
end

