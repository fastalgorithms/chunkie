function [eva_mats,sk_targ] = get_evamat_postproc_im(chnkr,clmparams,targs, ...
   targdomain,sk,eps,opts)
  %[~,ntarg] = size(targs);
  
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
    
    idomlist = 1:ndomain;
    if(isfield(opts,'idomlist'))
        idomlist = opts.idomlist;
    end
    
    chnkrtotal = merge(chnkr);
    r = chnkrtotal.r;
    npts = chnkrtotal.k*chnkrtotal.nch;
    r = reshape(r,[2,npts]);
    iiind = imag(r(1,:)) ~= 0;
    
    rnorms = normals(chnkrtotal);
    
       
    r_im = r(:,iiind);
    rnorms_im = rnorms(:,iiind);
    
    r_skel_im = r_im(:,sk);
    rnorms_skel_im = rnorms_im(:,sk);
    nmax = 2*numel(sk);
 
    
    eva_mats = cell(1,ndomain);
    sk_targ = cell(1,ndomain);
    for i=1:length(idomlist)
        idom = idomlist(i);
        if ~isempty(list{idom})
            
            srcinfo.r = r_skel_im;
            srcinfo.n = rnorms_skel_im;
            targinfo.r = targs(:,list{idom});
            ktmp = chnk.helm2d.kern(k(idom),srcinfo,targinfo,'evalg',coef(idom));
            bmat = [ktmp(:,:,1).';ktmp(:,:,2).';ktmp(:,:,3).';ktmp(:,:,4).';...
                ktmp(:,:,5).';ktmp(:,:,6).'];
            btmp2 = chnk.flam.rand_fft_transf(bmat,4*nmax);
            btmp3 = transpose(btmp2);
            btmp4 = chnk.flam.rand_fft_transf(btmp3,4*nmax);
            tol = eps;
            [SKT,~,~] = id(btmp4,tol); 

            bmat_red = btmp2(SKT(:),:);
            [SK_fin,RD_fin,T_fin] = id(bmat_red,tol); 

            eva_mat = zeros([numel(SK_fin),numel(targs(1,list{idom}))]);
            eva_mat(:,SK_fin) = eye(numel(SK_fin));
            eva_mat(:,RD_fin) = T_fin;
            eva_mat = transpose(eva_mat);
            eva_mats{idom} = eva_mat;
            sk_targ{idom} = SK_fin;
        end
    end
     
end