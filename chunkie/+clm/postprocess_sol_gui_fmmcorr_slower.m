function [u,gradu] = postprocess_sol_gui_fmmcorr_slower(chnkr,clmparams,targs,targdomain,eps0,sol,sk,exp_mat)
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
    srcinfo.charges = dens_c(irind).';
    
    r_im = r(:,iiind);
    rnorms_im = rnorms(:,iiind);
    
    r_skel_im = r_im(:,sk);
    rnorms_skel_im = rnorms_im(:,sk);
    dens_c_im = dens_c(iiind);
    dens_c_skel = exp_mat*dens_c_im;
    
    dens_d_im = dens_d(iiind);
    dens_d_skel = exp_mat*dens_d_im;
    
    pg = 0;
    pgt = 2;
    for i=1:ndomain
        if ~isempty(list{i})
            srcinfo.dipstr = dens_d(irind).'*coef(i);
            [UU] = hfmm2d(eps0,k(i),srcinfo,pg,targs(:,list{i}),pgt);
            u(list{i}) = UU.pottarg;
            gradu(:,list{i}) = UU.gradtarg;
            srcinfo.r = r_skel_im;
            srcinfo.n = rnorms_skel_im;
            targinfo.r = targs(:,list{i});
            ktmp = chnk.helm2d.kern(k(i),srcinfo,targinfo,'evalg',coef(i));
            %[~,ns] = size(srcinfo.r);
            u(list{i}) = u(list{i}) + (ktmp(:,:,1)*dens_d_skel).';
            u(list{i}) = u(list{i}) + (ktmp(:,:,2)*dens_c_skel).';
            gradu(1,list{i}) = gradu(1,list{i}) + (ktmp(:,:,3)*dens_d_skel).';
            gradu(1,list{i}) = gradu(1,list{i}) + (ktmp(:,:,4)*dens_c_skel).';
            gradu(2,list{i}) = gradu(2,list{i}) + (ktmp(:,:,5)*dens_d_skel).';
            gradu(2,list{i}) = gradu(2,list{i}) + (ktmp(:,:,6)*dens_c_skel).';
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
