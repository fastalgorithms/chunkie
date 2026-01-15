%%%%%%%%%%  SLOW CODE FOR FREDHOLM DETERMINANT USING CHEBYSHEV %%%%%%%%%
% m = 4;
% clear;
% f = @(zk) get_determinant(zk, get_chunker(1));
% |f - p| < 1.3/n^2
chnkreps = get_chunker(3);
f = chebfun(@(zk) get_determinant(zk, chnkreps), [3.3,3.5]);
plot(real(f), 'k.');
rts1 = roots(f, 'complex');


function det1 = get_determinant(zk, chnkr1)

    Dk = 2*kernel('helm', 'd', zk);  
    A = chunkermat(chnkr1, Dk);   
%%
    A = A + eye(chnkr1.npt);
%%
    det1 = det(A);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   CHEBYSHEV NODES  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ns = 1:5:200;
% err = zeros(size(ns));
% for jj = 1:length(ns)    
%     n = ns(jj);
%     t = 0:1:n;
%     x = cos(2*pi*t/(2*n + 1));              % Interpolating points
%     for kk = 1:length(t)
%         g(kk) = f(x(kk));                   % g = f(cos(t))
%     end
%     A = cos( (2 * pi / (2 * n + 1)) * [0:1:n].' * [0:1:n] );
%     d = 1/(2*n+1) * A * (g.*[1 2*ones(1,n)]).';
%     p = @(x) (cos(x*[0:1:n]).*[1 2*ones(1,n)])*d;
%     k = 10;
%     w = 0:pi/k:pi;
%     for ii = 1:length(w)
%         z(ii) = p(w(ii));
%         y(ii) = f(cos(w(ii)));
%     end
%     err(jj) = abs(max(y - z));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  FREDHOLM DETERMINANT  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% function det1 = get_determinant(zk, chnkr1)
% 
%     Dk = 2*kernel('helm', 'd', zk);  
%     A = chunkermat(chnkr1, Dk);   
% %%
%     A = A + eye(chnkr1.npt);
% %%
% 
%     det1 = det(A);
% end