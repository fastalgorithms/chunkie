function [varargout] = helm_axi_all(rs,drs,dzs,htables,ifun)

    r0s   = sqrt(rs.^2+(rs+drs).^2+dzs.^2);
    alphs = (drs.^2+dzs.^2)./r0s.^2;

    int   = zeros(size(alphs));

    iflag_rk = (ifun == 1 || ifun == 4 || ifun == 5);
    iflag_ik = (ifun == 2 || ifun == 4);
    iflag_dk = (ifun == 3 || ifun == 4 || ifun == 5);
    
    if (iflag_rk)
    rk = cell(6);
    for ii=1:6
        rk{ii} = int;
    end
    end

    if (iflag_ik)
    ik = cell(6);
    for ii=1:6
        ik{ii} = int;
    end
    end

    if (iflag_dk)
    dk = cell(6);
    for ii=1:6
        dk{ii} = int;
    end
    end
    
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
    
    
    if (iflag_rk)
    s ...
        = chnk.axissymhelm2d.helm_axi_close_table(r0s(iclose),alphs(iclose),1,htables);
    for ii=1:6
        rk{ii}(iclose) = s{ii};
    end
    end

    if (iflag_ik)
    s ...
        = chnk.axissymhelm2d.helm_axi_close_table(r0s(iclose),alphs(iclose),2,htables);
    for ii=1:6
        ik{ii}(iclose) = s{ii};
    end
    end

    if (iflag_dk)
    s ...
        = chnk.axissymhelm2d.helm_axi_close_table(r0s(iclose),alphs(iclose),3,htables);
    for ii=1:6
        dk{ii}(iclose) = s{ii};
    end
    end

    
 	s ...
        = chnk.axissymhelm2d.helm_axi_smooth(r0s(imidnearnear),alphs(imidnearnear),...
        ifun,htables.xlege_midnear,htables.wlege_midnear);
    
    if (iflag_rk)
    for ii=1:6
        rk{ii}(imidnearnear) = s.rk{ii};
    end
    end

    if(iflag_ik)
    for ii=1:6
        ik{ii}(imidnearnear) = s.ik{ii};
    end
    end

    if(iflag_dk)
    for ii=1:6
        dk{ii}(imidnearnear) = s.dk{ii};
    end
    end

	s ...
        = chnk.axissymhelm2d.helm_axi_smooth(r0s(imidnear),alphs(imidnear),...
        ifun,htables.xlege_midnear,htables.wlege_midnear);
    
    if (iflag_rk)
    for ii=1:6
        rk{ii}(imidnear) = s.rk{ii};
    end
    end

    if (iflag_ik)
    for ii=1:6
        ik{ii}(imidnear) = s.ik{ii};
    end
    end

    if (iflag_dk)
    for ii=1:6
        dk{ii}(imidnear) = s.dk{ii};
    end
    end

    s ...
        = chnk.axissymhelm2d.helm_axi_smooth(r0s(imid),alphs(imid),...
        ifun,htables.xlege_mid,htables.wlege_mid);

    if(iflag_rk)
    for ii=1:6
        rk{ii}(imid) = s.rk{ii};
    end
    end

    if (iflag_ik)
    for ii=1:6
        ik{ii}(imid) = s.ik{ii};
    end
    end

    if (iflag_dk)
    for ii=1:6
        dk{ii}(imid) = s.dk{ii};
    end
    end

    s ...
        = chnk.axissymhelm2d.helm_axi_smooth(r0s(ifar),alphs(ifar),...
        ifun,htables.xlege,htables.wlege);

    if (iflag_rk)
    for ii=1:6
        rk{ii}(ifar) = s.rk{ii};
    end
    end

    if (iflag_ik)
    for ii=1:6
        ik{ii}(ifar) = s.ik{ii};
    end
    end

    if (iflag_dk)
    for ii=1:6
        dk{ii}(ifar) = s.dk{ii};
    end
    end



    if(iflag_rk)
        [dsk] = proc_kerns(rs,drs,dzs,rk);
        varargout{1} = dsk;
        if (ifun == 4)
            [dsik] = proc_kerns(rs,drs,dzs,ik);
            [dsdiff] = proc_kerns(rs,drs,dzs,dk);
            varargout{2} = dsik;
            varargout{3} = dsdiff;
        end
    end
    if(ifun == 2)
        [dsik] = proc_kerns(rs,drs,dzs,ik);
        varargout{1} = dsik;
    end
    if(ifun==3)
        [dsdiff] = proc_kerns(rs,drs,dzs,dk);
        varargout{1} = dsdiff;
    end
    if(ifun==5)
        [dsdiff] = proc_kerns(rs,drs,dzs,dk);
        varargout{1} = dsdiff;
    end

end

function [ds] = proc_kerns(rs,drs,dzs,s)
    [dout] = chnk.axissymhelm2d.der_ak_to_grad(rs,drs,dzs,s{1},s{2},...
    s{3},s{6},s{5},s{4});
    [ds] = chnk.axissymhelm2d.div_by_kap(rs,drs,dzs,dout);
    ds.intdrz = -ds.intdrz;
     ds.intdzz = -ds.intdzz;
end