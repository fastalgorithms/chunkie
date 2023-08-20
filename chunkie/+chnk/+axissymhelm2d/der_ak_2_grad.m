function [dout] = der_ak_2_grad(r,q,z,i,ida,idk,...
    idaa,idak,idkk)

    dout = {};
    dout.int = i;
    
    ds = chnk.axissymhelm2d.dalph(r,q,z);
    
    intdr = ds.ar.*ida +ds.kr.*idk;
    intdq = ds.adr.*ida+ds.kdr.*idk;
    intdz = ds.az.*ida +ds.kz.*idk;
    
    intdrq= ds.ardr.*ida  +  ds.krdr.*idk  +  ...
            ds.ar.*ds.adr.*idaa +  ds.kr.*ds.kdr.*idkk +  ...
        ds.ar.*ds.kdr.*idak+ds.kr.*ds.adr.*idak;

    intdrz= ds.arz.*ida  +  ds.krz.*idk  +  ...
            ds.ar.*ds.az.*idaa +  ds.kr.*ds.kz.*idkk +  ...
        ds.ar.*ds.kz.*idak+ds.kr.*ds.az.*idak;    

    intdqz= ds.adrz.*ida  +  ds.kdrz.*idk  +  ...
            ds.adr.*ds.az.*idaa +  ds.kdr.*ds.kz.*idkk +  ...
        ds.adr.*ds.kz.*idak+ds.kdr.*ds.az.*idak;     
    
    intdzz= ds.azz.*ida  +  ds.kzz.*idk  +  ...
            ds.az.*ds.az.*idaa +  ds.kz.*ds.kz.*idkk +  ...
        ds.az.*ds.kz.*idak+ds.kz.*ds.az.*idak;    
    
    dout.intdr = intdr;
    dout.intdq = intdq;
    dout.intdz = intdz;
    
    dout.intdrq = intdrq;
    dout.intdrz = intdrz;
    dout.intdqz = intdqz;
    dout.intdzz = intdzz;
    
end

