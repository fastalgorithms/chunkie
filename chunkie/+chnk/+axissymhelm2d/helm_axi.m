function [ds] = helm_axi(r,dr,dz,ifun,htables)

    r0   = sqrt(r.^2+(r+dr).^2+dz.^2);
    alph = (dr.^2+dz.^2)./r0.^2;

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
    
    [dout] = chnk.axissymhelm2d.der_ak_2_grad(r,dr,dz,kern,kernda,kerndk,...
                 kerndaa,kerndak,kerndkk);
    [ds] = chnk.axissymhelm2d.div_by_kap(r,dr,dz,dout);
    ds.intdrz = -ds.intdrz;
    ds.intdzz = -ds.intdzz;

end

