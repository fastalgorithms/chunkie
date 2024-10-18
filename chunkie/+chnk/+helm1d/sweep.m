function pot = sweep(uin,inds,ts,wts,zkE)
%sweeping algorithm for preconditioner

    vexpsp = exp(1i*diff(ts)*zkE);

    u = zeros([numel(wts),1]);
    u(inds) = uin;

    nt = numel(u);
    chrg = u.*(wts); %got rid of scaling: %dm^2*u.*(wts)/zkE;
    
    voutp = zeros([nt,1]);
    voutp(1) = chrg(1);
    
    for i=2:nt
        voutp(i) = vexpsp(i-1)*voutp(i-1)+chrg(i);
    end
    
    voutm = zeros([nt,1]);
    chrg  = flipud(chrg);
    vexpsp= flipud(vexpsp);
    
    for i=2:nt
        voutm(i) = vexpsp(i-1)*(voutm(i-1)+chrg(i-1));
    end
    
    pot = voutp + flipud(voutm);
end