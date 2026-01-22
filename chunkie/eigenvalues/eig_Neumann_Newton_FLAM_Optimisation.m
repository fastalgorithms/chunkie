chnkr = get_chunker(4);
xx = chnkr.r(1,:);
yy = chnkr.r(2,:);

u = (xx.^2 - yy.^2 + yy).*sin(xx) + cos(xx).*(0.5*yy.^2 + 0.3*yy);
v = (xx - yy.^3 + xx.^2).*cos(-yy) + sin(yy).*(xx.^2 + xx);

% A = @(zk) get_matrix(zk, chnkr);
% derA = @(zk) get_dermatrix(zk, chnkr);
% y = @(zk) ones(1, length(u))*( ( (u.').*( A(zk)\(v.') ) ).*(chnkr.wts(:)) );
% dery = @(zk) ones(1, length(u))*...
%     (   (u.').*(  A(zk) \( derA(zk)*(A(zk)\(v.')) )  ).*(chnkr.wts(:))   );

F = @(zk) get_matrix(zk, chnkr);
derA = @(zk) get_dermatrix(zk, chnkr);
y = @(zk) ones(1, length(u))*( ( (u.').*( rskelf_sv(F(zk),v.') ) ).*(chnkr.wts(:)) );
dery = @(zk) ones(1, length(u))*...
    (   (u.').*(  rskelf_sv(F(zk), ( derA(zk)*(rskelf_sv(F(zk),v.')) ) )  ).*(chnkr.wts(:))   );

% rskelf_sv(F(zk),v)


fcm = @(zk) 1/y(zk);
fnr = @(zz) zz - y(zz)/dery(zz);    % defn of Newton-Raphson opti fn.

z0 = 3.409274696865717;       % x1 = 3.409274686865717; x2 = 3.409274686865787;
n = 10;% number of iterations
for i = 1:n
    z = fnr(z0);
    y2 = fcm(z);
    z0 = z;
    fprintf('i=%d  f=%d,  z0=%d\n',i, abs(y2), real(z0));
    if (abs(y2) < 1e-12)
        break;
    end
end
niters = i;

function [F] = get_matrix(zk, chnkr1)

    Dk = 2*kernel('helm', 'd', zk);  
    % A = chunkermat(chnkr1, Dk);
    F = chunkerflam(chnkr1, Dk, 1.0);
% %%
%     A = A + eye(chnkr1.npt);
end


function [derA] = get_dermatrix(zk, chnkr1)
    derDk = 2*kernel('helm', 'freq_diff', zk);  
    derA = chunkermat(chnkr1, derDk);   
end
