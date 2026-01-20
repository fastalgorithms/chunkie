%%%%%%%%     CENTRAL DIFFERENCE TEST FOR DERIVATIVE OF DETERMINANT   %%%%%%%%
clear;
chnkr = get_chunker(1);
det = @(zk) get_determinant(zk, chnkr);
derdet = @(zk) get_derdeterminant(zk, chnkr);

zk = 1:0.2:50;
err = zeros(1,length(zk));
err1 = zeros(1,length(zk));
h0 = 1e-4;
zz = 3.2;
for i = 1:length(zk)
     h = h0/zk(i);
     Dp = det(zz + h);  
     Dm = det(zz - h);
     derp = derdet(zz);
     err(i) = norm((Dp - Dm)/h - 2*derp);
end
loglog(zk, err, 'k.');
hold on; loglog(zk, 150*(h0./zk).^2, 'b.');


%%%%%%%%     NEWTON FOR DERIVATIVE OF DETERMINANT   %%%%%%%%





chnkreps = get_chunker(3.0);
det = @(zk) get_determinant(zk, chnkreps);
derdet = @(zk) get_derdeterminant(zk, chnkreps);
fnr = @(zz) zz - det(zz)/derdet(zz); % defn of Newton-Raphson opti fn.
%%
z0 = 3.4093; 
n = 10; % number of iterations
for i = 1:n
    [z, y] = update_iterate(z0, chnkreps);
    fprintf('i=%d  f=%d,  z0=%d\n',i, abs(y), real(z0));
    z0 = z;
    if (abs(y) < 1e-12)
        break;
    end
end
niters = i;


function [z, y] = update_iterate(zk, chnkr1)
    Dk = 2*kernel('helm', 'd', zk);  
    A = chunkermat(chnkr1, Dk, opts);   

    A = A + eye(chnkr1.npt);
    
    y = det(A);
    
    derDk = 2*kernel('helm', 'freq_diff', zk);  
    derA = chunkermat(chnkr1, derDk, opts);   
    z = zk - 1/trace(derA*inv(A));
end

function det1 = get_determinant(zk, chnkr1)

    Dk = 2*kernel('helm', 'd', zk);  
    A = chunkermat(chnkr1, Dk);   
%%
    A = A + eye(chnkr1.npt);
    
%%
    det1 = det(A);
    derdet1 = trace(det1*inv(A)*derA);
end


function derdet1 = get_derdeterminant(zk, chnkr1)

    Dk = 2*kernel('helm', 'd', zk);  
    A = chunkermat(chnkr1, Dk);   
    derDk = 2*kernel('helm', 'freq_diff', zk);  
    derA = chunkermat(chnkr1, derDk);   
%%
    A = A + eye(chnkr1.npt);
%%
    derdet1 = trace(det(A)*inv(A)*derA);
end


