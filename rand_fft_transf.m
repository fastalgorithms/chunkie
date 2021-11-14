function [atmp] = rand_fft_transf(amat,nout)
    sz = size(amat);
    nr = sz(1);
    nc = sz(2);
    dels = exp(1i*2*pi*rand([nr,1]));
    atmp = amat.*dels;
    %atmp = amat.*repmat(dels,[nr,1]);
    atmp = fft(atmp,nr)/sqrt(nr);
    is = randsample(nr,nout);
    atmp = atmp(is,:)*sqrt(nr/nout);
end

