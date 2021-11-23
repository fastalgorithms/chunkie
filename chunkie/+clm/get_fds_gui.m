function [Fskel,Fskel2,skel_struct,opts_perm,M,RG] = get_fds_gui(chnkr,clmparams,eps)
    nch = clmparams.nch;
    k = clmparams.ngl;
    np = sum(nch)*k;
    nn = 2*np;

    opts = [];
    opts.nonsmoothonly = true;

    [M,RG] = clm.get_mat_gui(chnkr,clmparams,opts);
    M = M + speye(nn);
    %tic, M1 = matfun(1:2*np,1:2*np); toc;
    

    
    opdims(1) = 2;
    opdims(2) = 2;
    
    chnkrtotal = merge(chnkr);
    

    rs = chnkrtotal.r(:,:);
    rs = repelem(rs,1,2);
    i_real = find((abs(imag(rs(1,:)))+abs(imag(rs(2,:)))) == 0);
    i_imag = find((abs(imag(rs(1,:)))+abs(imag(rs(2,:)))) > 0);
    iperm = [i_real,i_imag];

    opts_perm = [];
    opts_perm.iperm = iperm;
    invperm = 1:nn;
    invperm(iperm) = 1:nn;
    opts_perm.invperm = invperm;
    opts_perm.ns = [numel(i_real),numel(i_imag)];
    
    allt1 = @(s,t) chnk.helm2d.kern(1.0,s,t,'trans1',1);
    wts = weights(chnkrtotal);
    matfun = @(i,j) chnk.flam.kernbyindex(i,j,chnkrtotal,wts,allt1,opdims,M,...
        opts_perm);
    xflam = chnkrtotal.r(:,:);
    xflam = repelem(xflam,1,2);
    xflam = xflam(:,opts_perm.iperm);
    rank_or_tol = eps;
    occ = 400;
    pxyfun = [];
    opts = [];
    opts.lvlmax = 10;
    opts.verb = 1;
    tic, Fskel = rskelf(matfun,xflam(:,1:opts_perm.ns(1)),occ,rank_or_tol,pxyfun,opts); toc;

    opts_perm.n_offset = opts_perm.ns(1);
    matfun2 = @(i,j) chnk.flam.kernbyindex(i,j,chnkrtotal,wts,allt1,opdims,M,...
        opts_perm);
    irange = (opts_perm.ns(1)+1):(opts_perm.ns(1)+opts_perm.ns(2));
    xflam_i = [real(xflam(1,irange));imag(xflam(1,irange))];
    %xflam_i = xflam(:,irange);
    occ = 400;
    tic, Fskel2 = rskelf(matfun2,xflam_i,occ,rank_or_tol,pxyfun,opts); toc;
    
    opts_perm.n_offset = 0;

    %tic, M1 = matfun(1:nn,1:nn); toc;
    n1 = opts_perm.ns(1);
    n2 = opts_perm.ns(2);

    A12 = matfun(1:n1,(n1+1):nn);
    A21 = matfun((n1+1):nn,1:n1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nmax = 400;
    tol  = eps;
    [skel_struct] = chnk.flam.skel_2by2blk(Fskel,Fskel2,A12,A21,nmax,tol);


end