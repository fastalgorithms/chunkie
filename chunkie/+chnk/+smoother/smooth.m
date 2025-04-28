function [chnkr,varargout] = smooth(verts,opts)

    nv = size(verts,2);

    k = 16;
    if isfield(opts,'k')
       k = opts.k; 
    end
    
    kquad = 32;
    if isfield(opts,'kquad')
        kquad = 32;
    end
    
    nch  = 3;
    nchs = nch*ones(nv,1);

    if isfield(opts,'nchs')
        if (numel(opts.nchs) == 1)
            nchs = opts.nchs*ones(nv,1);
        elseif(size(opts.nchs,2) == nchs)
            nchs = opts.nchs;
        end
    end

    n_newton = 200;
    if isfield(opts,'n_newton')
        n_newton = opts.n_newton;
    end

    etol = 1E-6;
    if isfield(opts,'etol')
        etol = opts.etol;
    end

    lam = 2.5;
    if isfield(opts,'lam')
        lam = opts.lam;
    end

    umesh = chnk.smoother.get_umesh(verts);
    dmesh = chnk.smoother.get_mesh(umesh, nchs, k);
    qmesh = chnk.smoother.get_mesh(umesh, nchs, kquad);

    sig0 = sqrt(5)*max(umesh.lengths);
    if isfield(opts,'sig0')
        sig0 = opt.sig0;
    end

    dpx = dmesh.pseudo_normals(1,:).';
    dpy = dmesh.pseudo_normals(2,:).';

    opts_newt = [];
    opts_newt.scales = 1;
    opts_newt.levels = 0.5;

    opts_newt.step_fact = 1;
    if isfield(opts,'step_fact')
        opts_newt.step_fact = opts.step_fact;
    end

    [~, nd] = size(dmesh.r);
    h = zeros(nd,1);

    for i=1:n_newton
        [h,err,err_by_pt] = chnk.smoother.newt_step(h,umesh,dmesh,...
            qmesh,sig0,lam,opts_newt);
        rt = dmesh.r;
        rt(1,:) = rt(1,:) + (h.*dpx).';
        rt(2,:) = rt(2,:) + (h.*dpy).';
        if (err < etol)
            break
        end
    end

    if (err >= etol)
        error('ERROR: surface smoother failed to converge');
        return
    end

    chnkr = chunker.chunkerpoints(rt);
    if (nargout > 1)
        varargout{1} = err;
    end
    if (nargout > 2)
        varargout{2} = err_by_pt;
    end
end