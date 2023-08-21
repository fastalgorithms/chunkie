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

    
    for ii=1:numel(alphs)
    
    r  = rs(ii);
    dr = drs(ii);
    dz = dzs(ii);
    alph = alphs(ii);
    r0   = r0s(ii);
    
    if (alph <0.2)
        
        ik = ceil(r0/(pi));
        ia = ceil(-log(alph*5)/log(2));
        
        if (ia >111 || ik >152)
            % disp('Out of bounds');
            ds = [];
            ds.int = 0;
            ds.intdr = 0;
            ds.intdq = 0;
            ds.intdz = 0;
            ds.intdrq = 0;
            ds.intdrz = 0;
            ds.intdqz = 0;
            ds.intdzz = 0;
            return
        end
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
        
    else
        [kern,kernda,kerndk,kerndaa,kerndak,kerndkk] ...
            = chnk.axissymhelm2d.helm_axi_smooth(r0,alph,ifun,htables.xlege,htables.wlege);
    end
    
    kerns(ii)    = kern;
    kernsda(ii)  = kernda;
    kernsdk(ii)  = kerndk;
    kernsdak(ii) = kerndak;
    kernsdkk(ii) = kerndkk;
    kernsdaa(ii) = kerndaa;
    
    end
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

