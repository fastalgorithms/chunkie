function [sout] ...
        = helm_axi_close_table(r0s,alphs,ifun,htables)
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

        
        ik = ceil(r0/(pi));
        ia = ceil(-log(alph*5)/log(2));
        
        if (ia >111 || ik >152)
            kern    = 0;
            kernda  = 0;
            kerndk  = 0;
            kerndak = 0;
            kerndkk = 0;
            kerndaa = 0;
        else
        krel = (r0-pi*(ik-1))*2/pi-1;
        arel = ((alph - 0.2*2^(-ia))/(0.2*2^(-ia))-0.5)*2;
       
        tas = cos((0:(htables.ncheb-1))*acos(arel));
        tks = cos((0:(htables.ncheb-1))*acos(krel));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        vv = squeeze(htables.allvs(:,:,ifun,ia,ik));
        kern = tas*vv*tks.';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        vdk = squeeze(htables.alldk(:,:,ifun,ia,ik));
        kerndk = tas*vdk*tks.';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        vda = squeeze(htables.allda(:,:,ifun,ia,ik));
        kernda = tas*vda*tks.';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        vdak = squeeze(htables.alldak(:,:,ifun,ia,ik));
        kerndak = tas*vdak*tks.';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        vdkk = squeeze(htables.alldkk(:,:,ifun,ia,ik));
        kerndkk = tas*vdkk*tks.';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        vdaa = squeeze(htables.alldaa(:,:,ifun,ia,ik));
        kerndaa = tas*vdaa*tks.';
        
        end
        kerns(ii) = kern;
        kernsdk(ii) = kerndk;
        kernsda(ii) = kernda;
        kernsdaa(ii)= kerndaa;
        kernsdak(ii)= kerndak;
        kernsdkk(ii)= kerndkk;
        
    end    
    
    sout = {};
    sout{1} = kerns;
    sout{3} = kernsdk;
    sout{2} = kernsda;
    sout{4} = kernsdkk;
    sout{5} = kernsdak;
    sout{6} = kernsdaa;
    
end

