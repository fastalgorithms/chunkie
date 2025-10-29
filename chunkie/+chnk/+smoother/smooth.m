function [chnkr,varargout] = smooth(verts,opts)
%SMOOTH attempts to generate a smoothed curve from a polygon defined 
%   by a collection of vertices, using the method in:
%
%   "A fast boundary integral method for high-order 
%           multiscale mesh generation"
%
%   Felipe Vico, Leslie Greengard, Michael O'Neil, Manas Rachh
%
%   Though it should be noted that nothing about this code is 
%   expected to be fast.
%
% Syntax: chnkr = chnk.smoother.smooth(verts,opts)
%
% Input:
%   verts - a 2 by n matrix containing the x and y locations of the 
%           vertices of the original polygon. They are assumed to be 
%           positively orientated.
%   opts  - options structure. available options (default settings)
%           opts.k  - number of Legendre nodes per panel on the 
%                     output chunker object                     (16)
%           opts.kquad - number of quadrature nodes used by
%                     surface smoother in computing integrals. Must 
%                     be greater than k.                        (32)
%           opts.nchs  - the number of panels per edge in the final
%                     chunker. Can either by a single positive integer, 
%                     in which case each edge will have the same number 
%                     of panels, or a vector of length (# of vertices), 
%                     in which case the ith entry will specify the number 
%                     of panels on the edge joining the ith vertex to the
%                     (i+1)st vertex.                            (3)
%           opts.n_newton - the maximum number of steps of newton to 
%                     take when trying to smooth the curve.     (200)
%           opts.etol   - the error tolerance to be used when running 
%                     newton.                                   (1E-10)
%           opts.lam   - the smoothing parameter. Larger lambdas will 
%                     make the outputted curve less smooth. Small lambdas 
%                     increase the risk of the smoother failing to 
%                     converge. If the smoother fails, try running with
%                     a smaller lam.                            (10)
%           opts.step_fact - a parameter to multiply the step lengths 
%                     by when doing newton steps. If the smoother fails 
%                     to converge, then sometimes decreasing this can 
%                     help. Should always be <=1                (1)
%
% Output:
%   chnkr - the requested chunker containing the smoothed geometry
%   err (optional) - an approximate error of the curve smoother error
%   err_by_pt (optional) - a vector containing point by point errors
%
% The method works by convolving the polygon with a normalized Gaussian 
% (with a spatially dependent variance). The outputted curve ideally 
% is chosen to lie on the 1/2 level set of the resulting function. Errors
% here represent the difference between the values of this function at 
% the final points on the outputted curve and 1/2.
%

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
        elseif(numel(opts.nchs) == nv)
            nchs = opts.nchs;
        end
    end

    n_newton = 200;
    if isfield(opts,'n_newton')
        n_newton = opts.n_newton;
    end

    etol = 1E-10;
    if isfield(opts,'etol')
        etol = opts.etol;
    end

    lam = 10;
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

    rt = reshape(rt,2,k,[]);
    chnkr = chunkerpoints(rt);
    if (nargout > 1)
        varargout{1} = err;
    end
    if (nargout > 2)
        varargout{2} = err_by_pt;
    end
end
