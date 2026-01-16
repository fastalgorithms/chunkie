% % % % clear;
% % % % % ff = @(x) x - sin(x)/cos(x);
% % % % f = @(x) sin(x);
% % % % g = @(x) cos(x);
% % % % % g = @(x) 2*sign(x).*x;
% % % % x0 = 3.05;
% % % % n = 10;
% % % % for i = 1:n
% % % %     y = f(x0);
% % % %     dery = g(x0);
% % % %     % x = ff(x0);
% % % %     x = x0 - y/dery;
% % % %     x0 = x;
% % % %     if abs(sin(x))<1e-16
% % % %         break;
% % % %     end
% % % % end
% % % % x;
% % % % return


chnkr = get_chunker(3.35);
xx = chnkr.r(1,:);
yy = chnkr.r(2,:);
u = (xx.^2 - yy.^2 + yy).*sin(xx) + cos(xx).*(0.5*yy.^2 + 0.3*yy);
v = (xx - yy.^3 + xx.^2).*cos(-yy) + sin(yy).*(xx.^2 + xx);
A = @(zk) get_matrix(zk, chnkr);
derA = @(zk) get_dermatrix(zk, chnkr);
y = @(zk) ones(1, length(u))*( ( (u.').*( A(zk)\(v.') ) ).*(chnkr.wts(:)) );
dery = @(zk) ones(1, length(u))*...
    (   (u.').*(  A(zk) \( derA(zk)*(A(zk)\(v.')) )  ).*(chnkr.wts(:))   );
fcm = @(zk) 1/y(zk);
fnr = @(zz) zz - y(zz)/dery(zz); % defn of Newton-Raphson opti fn.

z0 = 3.409274696865717;% x1 = 3.409274686865717; x2 = 3.409274686865787;
n = 10;% number of iterations
for i = 1:n
    z = fnr(z0);
    y2 = fcm(z);
    z0 = z;
    if (abs(y2) < 1e-12)
        break;
    end
end
niters = i;


function [A] = get_matrix(zk, chnkr1)

    Dk = 2*kernel('helm', 'd', zk);  
    A = chunkermat(chnkr1, Dk);     
%%
    A = A + eye(chnkr1.npt);
end


function [derA] = get_dermatrix(zk, chnkr1)
    derDk = 2*kernel('helm', 'freq_diff', zk);  
    derA = chunkermat(chnkr1, derDk);   
end
