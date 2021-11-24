function [u,gradu] = postprocess_sol_gui_fmmcorr_slower_old(chnkr,clmparams,targs,targdomain,eps0,sol)
    [~,ntarg] = size(targs);
    u = zeros(1,ntarg);
    gradu = zeros(2,ntarg);
    

    if isfield(clmparams,'k')
      k = clmparams.k;
    end
    if isfield(clmparams,'ndomain')
      ndomain = clmparams.ndomain;
    end
    if isfield(clmparams,'coef')
        coef = clmparams.coef;
    end

    list = cell(1,ndomain);
    for i=1:ndomain
        list{i} = find(targdomain==i);
    end
    
    chnkrtotal = merge(chnkr);
    r = chnkrtotal.r;
    npts = chnkrtotal.k*chnkrtotal.nch;
    r = reshape(r,[2,npts]);
    
    irind = imag(r(1,:)) == 0;
    iiind = imag(r(1,:)) ~= 0;
    
    wts = weights(chnkrtotal);
    wts = wts(:);
    rnorms = normals(chnkrtotal);

    dens_d = sol(1:2:2*npts).*wts;
    dens_c = sol(2:2:2*npts).*wts;
    srcinfo.sources = r(:,irind);
    srcinfo.dipvec = rnorms(:,irind);
    srcinfo.charges = dens_c(irind);
    
    for i=1:ndomain
        if ~isempty(list{i})
            srcinfo.dipstr = dens_d(irind)*coef(i);
            [u(list{i}),gradu(:,list{i})] = hfmm2d(eps0,k(i),srcinfo,targs(:,list{i}));
            srcinfo.r = r(:,iiind);
            srcinfo.n = rnorms(:,iiind);
            targinfo.r = targs(:,list{i});
            ktmp = chnk.helm2d.kern(k(i),srcinfo,targinfo,'evalg',coef(i));
            %[~,ns] = size(srcinfo.r);
            u(list{i}) = u(list{i}) + (ktmp(:,:,1)*dens_d(iiind)).';
            u(list{i}) = u(list{i}) + (ktmp(:,:,2)*dens_c(iiind)).';
            gradu(1,list{i}) = gradu(1,list{i}) + (ktmp(:,:,3)*dens_d(iiind)).';
            gradu(1,list{i}) = gradu(1,list{i}) + (ktmp(:,:,4)*dens_c(iiind)).';
            gradu(2,list{i}) = gradu(2,list{i}) + (ktmp(:,:,5)*dens_d(iiind)).';
            gradu(2,list{i}) = gradu(2,list{i}) + (ktmp(:,:,6)*dens_c(iiind)).';
            j=i+1;
            if j > ndomain
              j = j - ndomain;
            end
        end
    end
    for i=1:ntarg
        if targdomain(i)==0
            u(i) = (u(i-1)+u(i+1))/2;
            gradu(:,i) = (gradu(:,i-1)+gradu(:,i+1))/2;
        end
    end
end
