function [SK,RD,exp_mat] = get_compressed_postproc_im(chnkr,clmparams)
    chnkrtotal = merge(chnkr);
    
    rs = chnkrtotal.r(:,:);
    %rs = repelem(rs,1,2);
    %i_real = find((abs(imag(rs(1,:)))+abs(imag(rs(2,:)))) == 0);
    i_imag = find((abs(imag(rs(1,:)))+abs(imag(rs(2,:)))) > 0);

    targs = clmparams.xylim;
    ks = clmparams.k;

    %ks = 1 + rand(100,1)*4;

    tol= 10^(-12);
    kmax = max(ks);
    %r_re = rs(:,i_real);
    r_im = rs(:,i_imag);

    ymin = targs(3);
    ymax = targs(4);
    xmax = targs(2);
    xmin = targs(1);

    x_used = [];
    y_used = [];


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    np = round((ymax-ymin)*kmax/(2*pi)*12);

    tys = (0:np)/np*(ymax-ymin) + ymin;
    txs = xmax*ones(size(tys));

    x_used = [x_used,txs];
    y_used = [y_used,tys];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tys = (0:np)/np*(ymax-ymin) + ymin;
    txs = xmin*ones(size(tys));

    x_used = [x_used,txs];
    y_used = [y_used,tys];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    np = round((xmax-xmin)*kmax/(2*pi)*12);

    txs = (0:np)/np*(xmax-xmin) + xmin;
    tys = ymin*ones(size(txs));

    x_used = [x_used,txs];
    y_used = [y_used,tys];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    np = 2*round((xmax-xmin)*kmax/(2*pi)*12);

    txs = (0:np)/np*(xmax-xmin) + xmin;
    tys = ymax*ones(size(txs));

    x_used = [x_used,txs];
    y_used = [y_used,tys];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nmax = 800;
    opts_ki2r =[];
    opts_ki2r.nmax  = nmax;
    opts_ki2r.ifcomp= true;

    xs = r_im(1,:);
    ys = r_im(2,:);
    xt = x_used;
    yt = y_used;

    atmp = [];

    for j=1:numel(ks)
        zk = ks(j);
        [atmp_j] = clm.get_kern_im2re(xs,ys,xt,yt,zk,opts_ki2r);
        atmp = [atmp;atmp_j];
        if (size(atmp,1)>size(atmp,2))
            sz = size(atmp,2);
            atmp = chnk.flam.rand_fft_transf(atmp,sz); 
        end    
    end

    btmp = chnk.flam.rand_fft_transf(atmp,nmax); 
    
    [SK,RD,T] = id(btmp,tol);

    exp_mat = zeros([numel(SK),numel(xs)]);
    exp_mat(:,SK) = eye(numel(SK));
    exp_mat(:,RD) = T;

end