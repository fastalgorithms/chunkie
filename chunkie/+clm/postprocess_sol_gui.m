function [u,gradu] = postprocess_sol_gui(chnkr,clmparams,targs,targdomain,eps0,sol)
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
    pg = 0;
    pgt = 2;
    for i=1:ndomain
        if ~isempty(list{i})
            r = chnkrtotal.r;
            npts = chnkrtotal.k*chnkrtotal.nch;
            r = reshape(r,[2,npts]);
            wts = weights(chnkrtotal);
            wts = wts(:);
            rnorms = normals(chnkrtotal);

            dens_d = sol(1:2:2*npts);
            dens_c = sol(2:2:2*npts);
            srcinfo.sources = r;
            srcinfo.dipvec = rnorms;
            srcinfo.charges = (dens_c.*wts).';
            srcinfo.dipstr = (dens_d.*wts*coef(i)).';
            UU = hfmm2d(eps0,k(i),srcinfo,pg,targs(:,list{i}),pgt);  
            u(list{i}) = UU.pottarg;
            gradu(:,list{i}) = UU.gradtarg;
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
