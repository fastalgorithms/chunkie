function [aout] = get_kern_im2re(xs,ys,xt,yt,zk,opts)
    if (nargin >5)
        if (isfield(opts,'nmax'))
            nmax = opts.nmax;
        else
            nmax = numel(xs);
        end
        if (isfield(opts,'ifcomp'))
            ifcomp = opts.ifcomp;
        else
            ifcomp = false;
        end
    end
    [Xt,Xs] = meshgrid(xt,xs);
    [Yt,Ys] = meshgrid(yt,ys);

    Z  = sqrt((Xt-Xs).^2+(Yt-Ys).^2);
    disp('Calculating Hankel functions')
    tic; [A1,A2] = hankm103(zk*Z); toc;
    %tic; A1 = besselh(0,1,zk*Z); toc;
    %tic; A2 = besselh(1,1,zk*Z); toc;
    
    A2a = A2./Z.*(Xt-Xs);
    A2b = A2./Z.*(Yt-Ys);
    A2 =[A2a,A2b];

    A12 = transpose([A1,A2]);
    if (nargin <6 || ifcomp==false)
        aout = A12;
        return;
    end
    size(A12)
    tic; aout = chnk.flam.rand_fft_transf(A12,nmax); toc;

end

