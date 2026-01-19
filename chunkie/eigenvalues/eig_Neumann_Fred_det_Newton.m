chnkreps = get_chunker(2.5);
det = @(zk) get_determinant(zk, chnkreps);
derdet = @(zk) get_derdeterminant(zk, chnkreps);
fnr = @(zz) zz - det(zz)/derdet(zz); % defn of Newton-Raphson opti fn.

z0 = 3; 
n = 5; % number of iterations
for i = 1:n
    z = fnr(z0);
    y = det(z);
    z0 = z;
    if (abs(y) < 1e-12)
        break;
    end
end
niters = i;


function det1 = get_determinant(zk, chnkr1)

    Dk = 2*kernel('helm', 'd', zk);  
    A = chunkermat(chnkr1, Dk);   
%%
    A = A + eye(chnkr1.npt);
    
%%
    det1 = det(A);
end


function derdet1 = get_derdeterminant(zk, chnkr1)

    Dk = 2*kernel('helm', 'd', zk);  
    A = chunkermat(chnkr1, Dk);   
    derDk = 2*kernel('helm', 'freq_diff', zk);  
    derA = chunkermat(chnkr1, derDk);   
%%
    A = A + eye(chnkr1.npt);
%%
    derdet1 = trace(adjoint(A)*derA);
end


